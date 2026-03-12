#!/usr/bin/env python3
"""
完整的工业级药物发现流程
Generation → ADMET → DrugLikeness → Docking → Ranking

流程:
1. 分子生成 (Diffusion Model)
2. 化学有效性验证
3. ADMET预测和筛选
4. 类药性评估和过滤
5. 分子对接评分
6. 综合排名和输出
"""

import torch
import numpy as np
import os
import json
import pandas as pd
from typing import List, Dict, Optional, Tuple
from pathlib import Path
from tqdm import tqdm
import logging
from dataclasses import dataclass, asdict

# 导入模块
from targetdiff_full import ModelConfig, EquivariantDiffusionModel, DiffusionSampler, load_model_checkpoint
from molecular_tools import PDBProcessor, MoleculeBuilder, MoleculeEvaluator, VinaDocking
from admet_filter import ADMETPredictor, DrugLikenessFilter, ADMETProfile, DrugLikenessProfile

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


@dataclass
class CompleteMoleculeProfile:
    """完整的分子评估档案"""
    # 基本信息
    mol_id: str
    smiles: str
    num_atoms: int
    num_bonds: int
    
    # 化学性质
    mw: float
    logp: float
    tpsa: float
    qed: float
    sa_score: float
    
    # ADMET
    admet_score: float
    caco2: float
    hia: float
    bbb: float
    ames_toxic: bool
    herg_risk: bool
    ld50: float
    
    # 类药性
    druglikeness_score: float
    lipinski_pass: bool
    lipinski_violations: int
    veber_pass: bool
    pains_alerts: int
    
    # 对接
    vina_score: float
    binding_affinity: float
    
    # 综合评分
    final_score: float
    rank: int
    recommended: bool


class IndustrialDrugDiscoveryPipeline:
    """工业级药物发现流程"""
    
    def __init__(self,
                 model_path: Optional[str] = None,
                 device: Optional[str] = None,
                 config: Optional[ModelConfig] = None):
        """初始化"""
        # 设备
        if device is None:
            self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        else:
            self.device = torch.device(device)
        
        logger.info(f"使用设备: {self.device}")
        
        # 配置
        self.config = config or ModelConfig()
        
        # 模型和工具
        self.model = load_model_checkpoint(model_path or 'pretrained_model.pt', self.config, self.device)
        self.sampler = DiffusionSampler(self.model, self.config, self.device)
        self.pdb_processor = PDBProcessor()
        self.mol_builder = MoleculeBuilder(self.config.atom_types)
        self.evaluator = MoleculeEvaluator()
        
        # 筛选器
        self.admet_predictor = ADMETPredictor()
        self.druglikeness_filter = DrugLikenessFilter()
        self.vina = VinaDocking()
        
        # 阈值设置
        self.thresholds = {
            'admet_score': 50.0,  # ADMET评分阈值
            'druglikeness_score': 50.0,  # 类药性评分阈值
            'lipinski_violations': 1,  # Lipinski违规数
            'pains_alerts': 0,  # PAINS警报数
            'vina_score': -6.0,  # Vina评分阈值(kcal/mol)
        }
        
        logger.info("工业级药物发现流程初始化完成")
    
    def run_complete_pipeline(self,
                             pdb_path: str,
                             num_candidates: int = 100,
                             num_atoms_range: tuple = (15, 35),
                             pocket_radius: float = 10.0,
                             sampling_steps: int = 500,
                             output_dir: str = 'industrial_results',
                             enable_docking: bool = True) -> List[CompleteMoleculeProfile]:
        """
        运行完整的工业级流程
        
        Args:
            pdb_path: PDB文件路径
            num_candidates: 初始候选分子数
            num_atoms_range: 原子数范围
            pocket_radius: 口袋半径
            sampling_steps: 采样步数
            output_dir: 输出目录
            enable_docking: 是否启用对接
        """
        os.makedirs(output_dir, exist_ok=True)
        
        logger.info("=" * 80)
        logger.info("工业级药物发现流程启动")
        logger.info("=" * 80)
        
        # 阶段1: 分子生成
        logger.info(f"\n阶段1/5: 生成 {num_candidates} 个候选分子")
        molecules = self._generate_molecules(
            pdb_path, num_candidates, num_atoms_range, 
            pocket_radius, sampling_steps, output_dir
        )
        logger.info(f"✓ 生成完成: {len(molecules)} 个有效分子")
        
        # 阶段2: ADMET预测
        logger.info(f"\n阶段2/5: ADMET预测和筛选")
        molecules = self._admet_screening(molecules)
        logger.info(f"✓ ADMET筛选后: {len(molecules)} 个分子")
        
        # 阶段3: 类药性评估
        logger.info(f"\n阶段3/5: 类药性评估和过滤")
        molecules = self._druglikeness_filtering(molecules)
        logger.info(f"✓ 类药性筛选后: {len(molecules)} 个分子")
        
        # 阶段4: 分子对接
        if enable_docking:
            logger.info(f"\n阶段4/5: 分子对接评分")
            molecules = self._molecular_docking(molecules, pdb_path, pocket_radius)
            logger.info(f"✓ 对接完成: {len(molecules)} 个分子")
        else:
            logger.info(f"\n阶段4/5: 跳过分子对接")
            for mol_data in molecules:
                mol_data['vina_score'] = 0.0
                mol_data['binding_affinity'] = 0.0
        
        # 阶段5: 综合评分和排名
        logger.info(f"\n阶段5/5: 综合评分和排名")
        profiles = self._calculate_final_scores(molecules)
        logger.info(f"✓ 评分完成: {len(profiles)} 个最终候选")
        
        # 保存结果
        self._save_results(profiles, output_dir)
        
        # 生成报告
        self._generate_report(profiles, output_dir)
        
        logger.info("\n" + "=" * 80)
        logger.info(f"完整流程完成! 结果保存在: {output_dir}")
        logger.info("=" * 80)
        
        return profiles
    
    def _generate_molecules(self, pdb_path, num_candidates, num_atoms_range, 
                           pocket_radius, sampling_steps, output_dir):
        """阶段1: 分子生成"""
        # 解析蛋白质
        protein_coords, _, _ = self.pdb_processor.parse_pdb(pdb_path)
        pocket_center = protein_coords.mean(axis=0)
        pocket_coords, _ = self.pdb_processor.extract_pocket(
            pdb_path, pocket_center, pocket_radius
        )
        
        # 转换为PyTorch
        protein_pos = torch.from_numpy(pocket_coords).float().to(self.device)
        protein_feat = self.model.atom_embedding(
            torch.zeros(len(pocket_coords), dtype=torch.long, device=self.device)
        )
        
        # 更新采样步数
        self.config.num_sampling_steps = sampling_steps
        
        molecules = []
        pbar = tqdm(range(num_candidates), desc="生成分子")
        
        for i in pbar:
            try:
                # 生成分子
                num_atoms = np.random.randint(num_atoms_range[0], num_atoms_range[1])
                pos, atom_type = self.sampler.sample(
                    num_atoms=num_atoms,
                    protein_pos=protein_pos,
                    protein_feat=protein_feat,
                    batch_size=1
                )
                
                # 构建RDKit分子
                pos = pos[0].cpu().numpy()
                atom_type = atom_type[0].cpu().numpy()
                mol = self.mol_builder.build_molecule(pos, atom_type, add_hydrogens=False)
                
                # 基本验证
                if mol is None or mol.GetNumAtoms() == 0:
                    continue
                
                # 保存
                mol_id = f"mol_{i:04d}"
                sdf_path = os.path.join(output_dir, 'molecules', f"{mol_id}.sdf")
                os.makedirs(os.path.dirname(sdf_path), exist_ok=True)
                
                try:
                    self.mol_builder.save_molecule(mol, sdf_path)
                except:
                    logger.warning(f"保存分子 {mol_id} 失败")
                    continue
                
                # 计算基本性质
                try:
                    smiles = Chem.MolToSmiles(mol)
                except:
                    smiles = "INVALID"
                
                molecules.append({
                    'mol_id': mol_id,
                    'mol': mol,
                    'smiles': smiles,
                    'pos': pos,
                    'atom_type': atom_type,
                    'sdf_path': sdf_path
                })
                
                pbar.set_postfix({'valid': len(molecules)})
                
            except Exception as e:
                logger.debug(f"生成分子 {i} 失败: {e}")
                continue
        
        return molecules
    
    def _admet_screening(self, molecules):
        """阶段2: ADMET预测"""
        filtered = []
        
        for mol_data in tqdm(molecules, desc="ADMET评估"):
            try:
                mol = mol_data['mol']
                
                # ADMET预测
                admet = self.admet_predictor.predict(mol)
                
                # 筛选标准
                if admet.overall_score < self.thresholds['admet_score']:
                    logger.debug(f"{mol_data['mol_id']}: ADMET评分过低 ({admet.overall_score:.1f})")
                    continue
                
                if admet.ames_toxicity:
                    logger.debug(f"{mol_data['mol_id']}: Ames致突变性阳性")
                    continue
                
                if admet.herg_inhibition:
                    logger.debug(f"{mol_data['mol_id']}: hERG抑制风险")
                    continue
                
                # 保存ADMET结果
                mol_data['admet'] = admet
                filtered.append(mol_data)
                
            except Exception as e:
                logger.debug(f"ADMET评估失败: {e}")
                continue
        
        return filtered
    
    def _druglikeness_filtering(self, molecules):
        """阶段3: 类药性评估"""
        filtered = []
        
        for mol_data in tqdm(molecules, desc="类药性评估"):
            try:
                mol = mol_data['mol']
                
                # 类药性评估
                druglikeness = self.druglikeness_filter.evaluate(mol)
                
                # 筛选标准
                if druglikeness.overall_score < self.thresholds['druglikeness_score']:
                    logger.debug(f"{mol_data['mol_id']}: 类药性评分过低")
                    continue
                
                if druglikeness.lipinski_violations > self.thresholds['lipinski_violations']:
                    logger.debug(f"{mol_data['mol_id']}: Lipinski违规过多")
                    continue
                
                if druglikeness.pains_alerts > self.thresholds['pains_alerts']:
                    logger.debug(f"{mol_data['mol_id']}: PAINS警报")
                    continue
                
                # 保存类药性结果
                mol_data['druglikeness'] = druglikeness
                filtered.append(mol_data)
                
            except Exception as e:
                logger.debug(f"类药性评估失败: {e}")
                continue
        
        return filtered
    
    def _molecular_docking(self, molecules, pdb_path, pocket_radius):
        """阶段4: 分子对接"""
        # 准备受体
        receptor_pdbqt = pdb_path.replace('.pdb', '.pdbqt')
        
        # 计算对接中心
        protein_coords, _, _ = self.pdb_processor.parse_pdb(pdb_path)
        pocket_center = protein_coords.mean(axis=0)
        
        for mol_data in tqdm(molecules, desc="分子对接"):
            try:
                mol = mol_data['mol']
                
                # 对接
                result = self.vina.dock(
                    receptor_pdbqt=receptor_pdbqt,
                    ligand_mol=mol,
                    center=tuple(pocket_center),
                    box_size=(pocket_radius * 2, pocket_radius * 2, pocket_radius * 2),
                    exhaustiveness=8
                )
                
                if result['success']:
                    mol_data['vina_score'] = result['affinity']
                    mol_data['binding_affinity'] = result['affinity']
                else:
                    mol_data['vina_score'] = 0.0
                    mol_data['binding_affinity'] = 0.0
                
            except Exception as e:
                logger.debug(f"对接失败: {e}")
                mol_data['vina_score'] = 0.0
                mol_data['binding_affinity'] = 0.0
        
        return molecules
    
    def _calculate_final_scores(self, molecules) -> List[CompleteMoleculeProfile]:
        """阶段5: 综合评分"""
        profiles = []
        
        for mol_data in molecules:
            try:
                mol = mol_data['mol']
                admet = mol_data.get('admet')
                druglikeness = mol_data.get('druglikeness')
                
                # 计算基本性质
                mw = Descriptors.MolWt(mol)
                logp = Descriptors.MolLogP(mol)
                tpsa = Descriptors.TPSA(mol)
                
                metrics = self.evaluator.calculate_metrics(mol)
                
                # 综合评分 (加权)
                final_score = (
                    admet.overall_score * 0.35 +  # ADMET 35%
                    druglikeness.overall_score * 0.25 +  # 类药性 25%
                    metrics.qed * 100 * 0.20 +  # QED 20%
                    min(abs(mol_data.get('vina_score', 0)) * 10, 100) * 0.20  # Vina 20%
                )
                
                # 推荐标准
                recommended = (
                    final_score >= 60 and
                    admet.overall_score >= 60 and
                    druglikeness.lipinski_pass and
                    mol_data.get('vina_score', 0) < -6.0
                )
                
                profile = CompleteMoleculeProfile(
                    mol_id=mol_data['mol_id'],
                    smiles=mol_data['smiles'],
                    num_atoms=mol.GetNumAtoms(),
                    num_bonds=mol.GetNumBonds(),
                    mw=mw,
                    logp=logp,
                    tpsa=tpsa,
                    qed=metrics.qed,
                    sa_score=metrics.sa_score,
                    admet_score=admet.overall_score,
                    caco2=admet.caco2_permeability,
                    hia=admet.hia,
                    bbb=admet.bbb_permeability,
                    ames_toxic=admet.ames_toxicity,
                    herg_risk=admet.herg_inhibition,
                    ld50=admet.ld50,
                    druglikeness_score=druglikeness.overall_score,
                    lipinski_pass=druglikeness.lipinski_pass,
                    lipinski_violations=druglikeness.lipinski_violations,
                    veber_pass=druglikeness.veber_pass,
                    pains_alerts=druglikeness.pains_alerts,
                    vina_score=mol_data.get('vina_score', 0.0),
                    binding_affinity=mol_data.get('binding_affinity', 0.0),
                    final_score=final_score,
                    rank=0,
                    recommended=recommended
                )
                
                profiles.append(profile)
                
            except Exception as e:
                logger.debug(f"评分失败: {e}")
                continue
        
        # 排名
        profiles.sort(key=lambda x: x.final_score, reverse=True)
        for i, profile in enumerate(profiles, 1):
            profile.rank = i
        
        return profiles
    
    def _save_results(self, profiles: List[CompleteMoleculeProfile], output_dir: str):
        """保存结果"""
        # JSON
        json_data = [asdict(p) for p in profiles]
        json_path = os.path.join(output_dir, 'complete_results.json')
        with open(json_path, 'w') as f:
            json.dump(json_data, f, indent=2)
        
        # CSV
        df = pd.DataFrame(json_data)
        csv_path = os.path.join(output_dir, 'complete_results.csv')
        df.to_csv(csv_path, index=False)
        
        # Top候选
        top_profiles = profiles[:20]
        top_df = pd.DataFrame([asdict(p) for p in top_profiles])
        top_csv = os.path.join(output_dir, 'top20_candidates.csv')
        top_df.to_csv(top_csv, index=False)
        
        logger.info(f"结果已保存: {output_dir}")
    
    def _generate_report(self, profiles: List[CompleteMoleculeProfile], output_dir: str):
        """生成分析报告"""
        report = {
            'summary': {
                'total_candidates': len(profiles),
                'recommended': sum(1 for p in profiles if p.recommended),
                'avg_final_score': np.mean([p.final_score for p in profiles]),
                'avg_admet_score': np.mean([p.admet_score for p in profiles]),
                'avg_druglikeness_score': np.mean([p.druglikeness_score for p in profiles]),
                'lipinski_pass_rate': sum(1 for p in profiles if p.lipinski_pass) / len(profiles) * 100,
                'avg_vina_score': np.mean([p.vina_score for p in profiles]),
            },
            'top_10': [asdict(p) for p in profiles[:10]],
            'statistics': {
                'mw_range': [min(p.mw for p in profiles), max(p.mw for p in profiles)],
                'logp_range': [min(p.logp for p in profiles), max(p.logp for p in profiles)],
                'qed_range': [min(p.qed for p in profiles), max(p.qed for p in profiles)],
            }
        }
        
        report_path = os.path.join(output_dir, 'analysis_report.json')
        with open(report_path, 'w') as f:
            json.dump(report, f, indent=2)
        
        # 打印摘要
        logger.info("\n" + "=" * 80)
        logger.info("分析摘要")
        logger.info("=" * 80)
        logger.info(f"总候选数: {report['summary']['total_candidates']}")
        logger.info(f"推荐数: {report['summary']['recommended']}")
        logger.info(f"平均综合评分: {report['summary']['avg_final_score']:.1f}")
        logger.info(f"平均ADMET评分: {report['summary']['avg_admet_score']:.1f}")
        logger.info(f"平均类药性评分: {report['summary']['avg_druglikeness_score']:.1f}")
        logger.info(f"Lipinski通过率: {report['summary']['lipinski_pass_rate']:.1f}%")
        logger.info(f"平均Vina评分: {report['summary']['avg_vina_score']:.2f} kcal/mol")


def main():
    """主函数"""
    import argparse
    
    parser = argparse.ArgumentParser(description='工业级药物发现流程')
    parser.add_argument('--pdb', type=str, required=True, help='PDB文件')
    parser.add_argument('--candidates', type=int, default=100, help='候选分子数')
    parser.add_argument('--num_atoms', type=int, nargs=2, default=[15, 35], help='原子数范围')
    parser.add_argument('--sampling_steps', type=int, default=500, help='采样步数')
    parser.add_argument('--no_docking', action='store_true', help='禁用对接')
    parser.add_argument('--output', type=str, default='industrial_results', help='输出目录')
    parser.add_argument('--device', type=str, default=None, help='设备')
    
    args = parser.parse_args()
    
    # 运行
    pipeline = IndustrialDrugDiscoveryPipeline(device=args.device)
    
    profiles = pipeline.run_complete_pipeline(
        pdb_path=args.pdb,
        num_candidates=args.candidates,
        num_atoms_range=tuple(args.num_atoms),
        sampling_steps=args.sampling_steps,
        output_dir=args.output,
        enable_docking=not args.no_docking
    )
    
    print(f"\n✓ 完成! 发现 {sum(1 for p in profiles if p.recommended)} 个推荐候选分子")


if __name__ == '__main__':
    main()
