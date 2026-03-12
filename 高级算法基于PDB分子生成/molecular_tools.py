#!/usr/bin/env python3
"""
分子工具模块
包含PDB处理、分子构建、Vina对接等功能
"""

import torch
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Crippen, Lipinski, QED
from rdkit.Chem import rdMolTransforms
from typing import List, Tuple, Optional, Dict
import os
import subprocess
import tempfile
import logging
from dataclasses import dataclass
import requests
from io import StringIO

logger = logging.getLogger(__name__)


@dataclass
class MoleculeMetrics:
    """分子评估指标"""
    qed: float  # 药物相似性
    sa_score: float  # 合成可及性
    logp: float  # 脂溶性
    mw: float  # 分子量
    hbd: int  # 氢键供体
    hba: int  # 氢键受体
    tpsa: float  # 极性表面积
    rotatable_bonds: int  # 可旋转键
    validity: bool  # 是否有效
    uniqueness: bool = True
    diversity: float = 1.0


class PDBProcessor:
    """PDB文件处理器"""
    
    @staticmethod
    def download_pdb(pdb_id: str, output_path: str) -> bool:
        """
        从RCSB PDB下载蛋白质结构
        Args:
            pdb_id: PDB ID (例如: '1h36')
            output_path: 输出文件路径
        """
        try:
            url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
            response = requests.get(url, timeout=30)
            
            if response.status_code == 200:
                with open(output_path, 'w') as f:
                    f.write(response.text)
                logger.info(f"成功下载PDB: {pdb_id} -> {output_path}")
                return True
            else:
                logger.error(f"下载失败: {pdb_id}, 状态码: {response.status_code}")
                return False
        except Exception as e:
            logger.error(f"下载PDB出错: {e}")
            return False
    
    @staticmethod
    def parse_pdb(pdb_path: str, extract_protein: bool = True) -> Tuple[np.ndarray, List[str], List[str]]:
        """
        解析PDB文件
        Args:
            pdb_path: PDB文件路径
            extract_protein: 是否只提取蛋白质原子
        Returns:
            coords: [N, 3] 坐标
            atoms: [N] 原子类型
            residues: [N] 残基信息
        """
        coords = []
        atoms = []
        residues = []
        
        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith('ATOM') or (not extract_protein and line.startswith('HETATM')):
                    try:
                        atom_name = line[12:16].strip()
                        res_name = line[17:20].strip()
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        element = line[76:78].strip()
                        
                        if not element:
                            # 从原子名推断元素
                            element = atom_name[0]
                        
                        coords.append([x, y, z])
                        atoms.append(element)
                        residues.append(res_name)
                    except (ValueError, IndexError):
                        continue
        
        if len(coords) == 0:
            raise ValueError(f"无法从 {pdb_path} 读取有效坐标")
        
        return np.array(coords), atoms, residues
    
    @staticmethod
    def extract_pocket(pdb_path: str, ligand_center: np.ndarray, pocket_radius: float = 10.0) -> Tuple[np.ndarray, List[str]]:
        """
        提取蛋白质口袋
        Args:
            pdb_path: PDB文件路径
            ligand_center: 配体中心坐标
            pocket_radius: 口袋半径(Å)
        Returns:
            pocket_coords: [N, 3] 口袋坐标
            pocket_atoms: [N] 原子类型
        """
        coords, atoms, residues = PDBProcessor.parse_pdb(pdb_path)
        
        # 计算到配体中心的距离
        distances = np.linalg.norm(coords - ligand_center, axis=1)
        
        # 筛选口袋内的原子
        mask = distances <= pocket_radius
        pocket_coords = coords[mask]
        pocket_atoms = [atoms[i] for i in range(len(atoms)) if mask[i]]
        
        logger.info(f"提取口袋: {len(pocket_coords)} 个原子 (半径 {pocket_radius}Å)")
        return pocket_coords, pocket_atoms
    
    @staticmethod
    def save_pdb(coords: np.ndarray, atoms: List[str], output_path: str, molecule_name: str = "LIG"):
        """保存为PDB文件"""
        with open(output_path, 'w') as f:
            f.write(f"HEADER    {molecule_name}\n")
            for i, (coord, atom) in enumerate(zip(coords, atoms), 1):
                f.write(f"HETATM{i:5d}  {atom:<3s} {molecule_name} A   1    "
                       f"{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}  1.00  0.00          {atom:>2s}\n")
            f.write("END\n")


class MoleculeBuilder:
    """分子构建器"""
    
    def __init__(self, atom_types: List[str]):
        self.atom_types = atom_types
        self.bond_length_dict = {
            ('C', 'C'): 1.54, ('C', 'N'): 1.47, ('C', 'O'): 1.43,
            ('C', 'S'): 1.82, ('C', 'F'): 1.35, ('C', 'P'): 1.85,
            ('C', 'Cl'): 1.77, ('C', 'Br'): 1.94, ('C', 'I'): 2.14,
            ('N', 'N'): 1.45, ('N', 'O'): 1.40, ('N', 'S'): 1.68,
            ('O', 'O'): 1.48, ('O', 'S'): 1.58, ('O', 'P'): 1.63,
            ('S', 'S'): 2.05
        }
    
    def build_molecule(self, coords: np.ndarray, atom_types: np.ndarray, add_hydrogens: bool = False) -> Chem.Mol:
        """
        从坐标和原子类型构建RDKit分子
        Args:
            coords: [N, 3] 原子坐标
            atom_types: [N] 原子类型索引
            add_hydrogens: 是否添加氢原子
        """
        mol = Chem.RWMol()
        
        # 添加原子
        atom_map = {}
        for i, atom_idx in enumerate(atom_types):
            if isinstance(atom_idx, torch.Tensor):
                atom_idx = atom_idx.item()
            
            if atom_idx >= len(self.atom_types):
                atom_idx = 0  # 默认为碳
            
            atom_symbol = self.atom_types[atom_idx]
            atom = Chem.Atom(atom_symbol)
            idx = mol.AddAtom(atom)
            atom_map[i] = idx
        
        # 添加构象
        conf = Chem.Conformer(len(atom_types))
        for i, coord in enumerate(coords):
            if isinstance(coord, torch.Tensor):
                coord = coord.cpu().numpy()
            conf.SetAtomPosition(i, tuple(coord.astype(float)))
        mol.AddConformer(conf)
        
        # 基于距离添加键
        self._add_bonds_by_distance(mol, coords)
        
        # 转换为不可变分子
        mol = mol.GetMol()
        
        # 清理分子 - 使用更宽松的选项
        try:
            # 尝试标准清理
            Chem.SanitizeMol(mol, catchErrors=True)
        except Exception as e:
            logger.warning(f"标准清理失败: {e}, 尝试宽松清理")
            try:
                # 宽松清理 - 跳过某些检查
                Chem.SanitizeMol(mol, 
                                sanitizeOps=Chem.SanitizeFlags.SANITIZE_FINDRADICALS|
                                           Chem.SanitizeFlags.SANITIZE_KEKULIZE|
                                           Chem.SanitizeFlags.SANITIZE_SETAROMATICITY|
                                           Chem.SanitizeFlags.SANITIZE_SETCONJUGATION|
                                           Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION,
                                catchErrors=True)
            except:
                logger.warning("宽松清理也失败，使用原始分子")
        
        # 尝试添加氢原子（通常不需要）
        if add_hydrogens:
            try:
                mol = Chem.AddHs(mol, addCoords=True)
            except:
                logger.warning("添加氢原子失败")
        
        return mol
    
    def _add_bonds_by_distance(self, mol: Chem.RWMol, coords: np.ndarray):
        """基于距离添加键"""
        n_atoms = len(coords)
        
        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                dist = np.linalg.norm(coords[i] - coords[j])
                
                atom_i = mol.GetAtomWithIdx(i).GetSymbol()
                atom_j = mol.GetAtomWithIdx(j).GetSymbol()
                
                # 获取标准键长
                bond_key = tuple(sorted([atom_i, atom_j]))
                expected_dist = self.bond_length_dict.get(bond_key, 1.5)
                
                # 判断是否成键
                if dist < expected_dist * 1.2:  # 20%容差
                    try:
                        if dist < expected_dist * 0.9:
                            mol.AddBond(i, j, Chem.BondType.DOUBLE)
                        else:
                            mol.AddBond(i, j, Chem.BondType.SINGLE)
                    except:
                        pass
    
    def save_molecule(self, mol: Chem.Mol, output_path: str):
        """保存分子"""
        ext = os.path.splitext(output_path)[1].lower()
        
        if ext == '.sdf':
            writer = Chem.SDWriter(output_path)
            writer.write(mol)
            writer.close()
        elif ext == '.mol2':
            Chem.MolToMolFile(mol, output_path.replace('.mol2', '.mol'))
            logger.info(f"分子已保存为MOL格式: {output_path.replace('.mol2', '.mol')}")
        elif ext == '.pdb':
            Chem.MolToPDBFile(mol, output_path)
        else:
            raise ValueError(f"不支持的文件格式: {ext}")
        
        logger.info(f"分子已保存: {output_path}")


class MoleculeEvaluator:
    """分子评估器"""
    
    @staticmethod
    def calculate_metrics(mol: Chem.Mol) -> MoleculeMetrics:
        """计算分子指标"""
        try:
            # 先尝试更新分子属性
            try:
                Chem.SanitizeMol(mol, catchErrors=True)
            except:
                pass
            
            # QED (药物相似性)
            try:
                qed_score = QED.qed(mol)
            except:
                qed_score = 0.0
            
            # 合成可及性 (近似)
            sa_score = MoleculeEvaluator._calculate_sa_score(mol)
            
            # 理化性质 - 使用try-except保护每个计算
            try:
                logp = Crippen.MolLogP(mol)
            except:
                logp = 0.0
                
            try:
                mw = Descriptors.MolWt(mol)
            except:
                # 手动计算分子量
                mw = sum([atom.GetMass() for atom in mol.GetAtoms()])
                
            try:
                hbd = Lipinski.NumHDonors(mol)
            except:
                hbd = 0
                
            try:
                hba = Lipinski.NumHAcceptors(mol)
            except:
                hba = 0
                
            try:
                tpsa = Descriptors.TPSA(mol)
            except:
                tpsa = 0.0
                
            try:
                rotatable_bonds = Lipinski.NumRotatableBonds(mol)
            except:
                rotatable_bonds = 0
            
            validity = True
        except Exception as e:
            logger.warning(f"计算分子指标失败: {e}")
            qed_score = 0.0
            sa_score = 10.0
            logp = mw = hbd = hba = tpsa = rotatable_bonds = 0
            validity = False
        
        return MoleculeMetrics(
            qed=qed_score,
            sa_score=sa_score,
            logp=logp,
            mw=mw,
            hbd=hbd,
            hba=hba,
            tpsa=tpsa,
            rotatable_bonds=rotatable_bonds,
            validity=validity
        )
    
    @staticmethod
    def _calculate_sa_score(mol: Chem.Mol) -> float:
        """
        计算合成可及性分数 (SA Score)
        简化版本,范围1-10,越小越容易合成
        """
        # 复杂度因素
        n_atoms = mol.GetNumAtoms()
        n_chiral = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
        n_rings = Chem.GetSSSR(mol)
        n_heteroatoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in [1, 6])
        
        # 简化计算
        complexity = (
            0.1 * n_atoms +
            0.5 * n_chiral +
            0.3 * n_rings +
            0.2 * n_heteroatoms
        )
        
        # 归一化到1-10
        sa_score = min(10, max(1, 1 + complexity / 10))
        return sa_score
    
    @staticmethod
    def check_lipinski(mol: Chem.Mol) -> Dict[str, bool]:
        """Lipinski五规则检查"""
        try:
            # 先尝试清理分子
            try:
                Chem.SanitizeMol(mol, catchErrors=True)
            except:
                pass
            
            # 安全计算每个属性
            try:
                mw = Descriptors.MolWt(mol)
            except:
                mw = sum([atom.GetMass() for atom in mol.GetAtoms()])
            
            try:
                logp = Crippen.MolLogP(mol)
            except:
                logp = 0.0
            
            try:
                hbd = Lipinski.NumHDonors(mol)
            except:
                hbd = 0
            
            try:
                hba = Lipinski.NumHAcceptors(mol)
            except:
                hba = 0
            
            return {
                'mw': mw <= 500,
                'logp': logp <= 5,
                'hbd': hbd <= 5,
                'hba': hba <= 10,
                'pass': all([mw <= 500, logp <= 5, hbd <= 5, hba <= 10])
            }
        except Exception as e:
            logger.warning(f"Lipinski检查失败: {e}")
            return {
                'mw': False,
                'logp': False,
                'hbd': False,
                'hba': False,
                'pass': False
            }


class VinaDocking:
    """
    AutoDock Vina对接
    需要安装: pip install vina meeko
    """
    
    def __init__(self, vina_executable: str = 'vina'):
        self.vina_executable = vina_executable
        self._check_vina()
    
    def _check_vina(self):
        """检查Vina是否可用"""
        try:
            result = subprocess.run([self.vina_executable, '--version'], 
                                   capture_output=True, text=True, timeout=5)
            if result.returncode == 0:
                logger.info(f"Vina可用: {result.stdout.strip()}")
            else:
                logger.warning("Vina不可用,对接功能将被禁用")
        except:
            logger.warning("Vina不可用,对接功能将被禁用")
    
    def prepare_receptor(self, pdb_path: str, output_pdbqt: str):
        """准备受体(需要MGLTools)"""
        try:
            cmd = f"prepare_receptor4.py -r {pdb_path} -o {output_pdbqt}"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=60)
            if result.returncode == 0:
                logger.info(f"受体准备完成: {output_pdbqt}")
                return True
            else:
                logger.error(f"受体准备失败: {result.stderr}")
                return False
        except Exception as e:
            logger.error(f"准备受体出错: {e}")
            return False
    
    def dock(self, 
             receptor_pdbqt: str,
             ligand_mol: Chem.Mol,
             center: Tuple[float, float, float],
             box_size: Tuple[float, float, float] = (20, 20, 20),
             exhaustiveness: int = 8,
             num_modes: int = 9) -> Dict:
        """
        执行对接
        Args:
            receptor_pdbqt: 受体PDBQT文件
            ligand_mol: 配体分子
            center: 对接盒子中心 (x, y, z)
            box_size: 对接盒子大小 (size_x, size_y, size_z)
            exhaustiveness: 搜索精度
            num_modes: 生成构象数
        Returns:
            结果字典,包含affinity等
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            # 保存配体
            ligand_pdbqt = os.path.join(tmpdir, 'ligand.pdbqt')
            self._mol_to_pdbqt(ligand_mol, ligand_pdbqt)
            
            # 对接输出
            output_pdbqt = os.path.join(tmpdir, 'output.pdbqt')
            
            # 构建Vina命令
            cmd = [
                self.vina_executable,
                '--receptor', receptor_pdbqt,
                '--ligand', ligand_pdbqt,
                '--out', output_pdbqt,
                '--center_x', str(center[0]),
                '--center_y', str(center[1]),
                '--center_z', str(center[2]),
                '--size_x', str(box_size[0]),
                '--size_y', str(box_size[1]),
                '--size_z', str(box_size[2]),
                '--exhaustiveness', str(exhaustiveness),
                '--num_modes', str(num_modes)
            ]
            
            try:
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
                
                if result.returncode == 0:
                    # 解析结果
                    affinity = self._parse_vina_output(result.stdout)
                    return {
                        'success': True,
                        'affinity': affinity,
                        'output': result.stdout
                    }
                else:
                    return {
                        'success': False,
                        'error': result.stderr
                    }
            except Exception as e:
                logger.error(f"对接出错: {e}")
                return {'success': False, 'error': str(e)}
    
    def _mol_to_pdbqt(self, mol: Chem.Mol, output_path: str):
        """将RDKit分子转为PDBQT"""
        try:
            from meeko import MoleculePreparation
            
            preparator = MoleculePreparation()
            preparator.prepare(mol)
            preparator.write_pdbqt_file(output_path)
        except Exception as e:
            logger.error(f"转换PDBQT失败: {e}")
            # 降级方案:转为PDB
            pdb_path = output_path.replace('.pdbqt', '.pdb')
            Chem.MolToPDBFile(mol, pdb_path)
    
    def _parse_vina_output(self, output: str) -> float:
        """解析Vina输出获取亲和力"""
        for line in output.split('\n'):
            if 'REMARK VINA RESULT' in line or line.strip().startswith('1'):
                parts = line.split()
                try:
                    return float(parts[1])
                except:
                    continue
        return 0.0


if __name__ == '__main__':
    # 测试代码
    logging.basicConfig(level=logging.INFO)
    
    # 测试PDB下载
    # PDBProcessor.download_pdb('1h36', 'test_1h36.pdb')
    
    # 测试分子构建
    coords = np.random.randn(10, 3) * 2
    atom_types = np.random.randint(0, 3, 10)
    
    builder = MoleculeBuilder(['C', 'N', 'O'])
    mol = builder.build_molecule(coords, atom_types)
    
    print(f"构建分子: {Chem.MolToSmiles(mol)}")
    print(f"原子数: {mol.GetNumAtoms()}")
    
    # 测试评估
    evaluator = MoleculeEvaluator()
    metrics = evaluator.calculate_metrics(mol)
    print(f"QED: {metrics.qed:.3f}, MW: {metrics.mw:.1f}")