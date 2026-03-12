#!/usr/bin/env python3
"""
ADMET预测和类药性筛选模块
实现完整的药物筛选流程

ADMET评估:
- Absorption: 吸收性
- Distribution: 分布性
- Metabolism: 代谢性
- Excretion: 排泄性
- Toxicity: 毒性

类药性规则:
- Lipinski规则 (Ro5)
- Veber规则
- Ghose规则
- Egan规则
- PAINS过滤
"""

import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, AllChem, QED
from rdkit.Chem import rdMolDescriptors, Fragments, MolSurf
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
import logging
import warnings

warnings.filterwarnings('ignore')
logger = logging.getLogger(__name__)


@dataclass
class ADMETProfile:
    """ADMET评估结果"""
    # Absorption (吸收)
    caco2_permeability: float  # Caco-2通透性 (-8到-4: 差, -4到-5.5: 中等, >-5.5: 好)
    hia: float  # 人肠道吸收 (0-100%)
    pgp_substrate: bool  # P-糖蛋白底物
    
    # Distribution (分布)
    vd: float  # 分布体积 (L/kg)
    bbb_permeability: float  # 血脑屏障通透性
    ppb: float  # 血浆蛋白结合率 (%)
    
    # Metabolism (代谢)
    cyp_substrate: Dict[str, bool]  # CYP酶底物
    cyp_inhibitor: Dict[str, bool]  # CYP酶抑制剂
    
    # Excretion (排泄)
    clearance: float  # 清除率 (mL/min/kg)
    half_life: float  # 半衰期 (小时)
    
    # Toxicity (毒性)
    ames_toxicity: bool  # Ames致突变性
    herg_inhibition: bool  # hERG抑制
    hepatotoxicity: bool  # 肝毒性
    skin_sensitization: bool  # 皮肤敏感
    ld50: float  # 半数致死量 (mg/kg)
    
    # 综合评分
    overall_score: float  # 0-100


@dataclass
class DrugLikenessProfile:
    """类药性评估结果"""
    # Lipinski规则 (Ro5)
    lipinski_pass: bool
    lipinski_violations: int
    
    # Veber规则
    veber_pass: bool
    
    # Ghose规则
    ghose_pass: bool
    
    # Egan规则
    egan_pass: bool
    
    # 其他规则
    lead_like: bool  # 先导化合物相似性
    fragment_like: bool  # 片段相似性
    
    # PAINS
    pains_alerts: int  # PAINS警报数
    
    # 综合评分
    overall_score: float  # 0-100


class ADMETPredictor:
    """
    ADMET预测器
    使用基于分子描述符的经验模型
    """
    
    def __init__(self):
        self.cyp_families = ['1A2', '2C9', '2C19', '2D6', '3A4']
    
    def predict(self, mol: Chem.Mol) -> ADMETProfile:
        """完整ADMET预测"""
        try:
            # 计算分子描述符
            descriptors = self._calculate_descriptors(mol)
            
            # Absorption
            caco2 = self._predict_caco2(descriptors)
            hia = self._predict_hia(descriptors)
            pgp = self._predict_pgp_substrate(descriptors)
            
            # Distribution
            vd = self._predict_vd(descriptors)
            bbb = self._predict_bbb(descriptors)
            ppb = self._predict_ppb(descriptors)
            
            # Metabolism
            cyp_sub = self._predict_cyp_substrate(descriptors)
            cyp_inh = self._predict_cyp_inhibitor(descriptors)
            
            # Excretion
            clearance = self._predict_clearance(descriptors)
            half_life = self._predict_half_life(descriptors)
            
            # Toxicity
            ames = self._predict_ames(descriptors)
            herg = self._predict_herg(descriptors)
            hepato = self._predict_hepatotoxicity(descriptors)
            skin = self._predict_skin_sensitization(descriptors)
            ld50 = self._predict_ld50(descriptors)
            
            # 综合评分
            overall = self._calculate_overall_score(
                caco2, hia, bbb, ames, herg, hepato, ld50
            )
            
            return ADMETProfile(
                caco2_permeability=caco2,
                hia=hia,
                pgp_substrate=pgp,
                vd=vd,
                bbb_permeability=bbb,
                ppb=ppb,
                cyp_substrate=cyp_sub,
                cyp_inhibitor=cyp_inh,
                clearance=clearance,
                half_life=half_life,
                ames_toxicity=ames,
                herg_inhibition=herg,
                hepatotoxicity=hepato,
                skin_sensitization=skin,
                ld50=ld50,
                overall_score=overall
            )
            
        except Exception as e:
            logger.error(f"ADMET预测失败: {e}")
            return self._get_default_profile()
    
    def _calculate_descriptors(self, mol: Chem.Mol) -> Dict:
        """计算分子描述符"""
        desc = {}
        
        try:
            # 基本性质
            desc['mw'] = Descriptors.MolWt(mol)
            desc['logp'] = Crippen.MolLogP(mol)
            desc['tpsa'] = Descriptors.TPSA(mol)
            desc['hbd'] = Lipinski.NumHDonors(mol)
            desc['hba'] = Lipinski.NumHAcceptors(mol)
            desc['rotatable_bonds'] = Lipinski.NumRotatableBonds(mol)
            desc['aromatic_rings'] = Descriptors.NumAromaticRings(mol)
            
            # 拓扑描述符
            desc['heavy_atoms'] = Lipinski.HeavyAtomCount(mol)
            desc['rings'] = Descriptors.RingCount(mol)
            desc['sp3_fraction'] = Descriptors.FractionCsp3(mol)
            
            # 表面积
            try:
                desc['labuteasa'] = Descriptors.LabuteASA(mol)
            except:
                desc['labuteasa'] = 0
            
            # 极性
            desc['num_hdonors'] = desc['hbd']
            desc['num_hacceptors'] = desc['hba']
            
        except Exception as e:
            logger.warning(f"描述符计算部分失败: {e}")
        
        return desc
    
    def _predict_caco2(self, desc: Dict) -> float:
        """
        预测Caco-2通透性
        经验公式基于LogP和TPSA
        """
        logp = desc.get('logp', 0)
        tpsa = desc.get('tpsa', 0)
        
        # 经验公式
        caco2 = -5.0 + 0.5 * logp - 0.01 * tpsa
        return np.clip(caco2, -8, -3)
    
    def _predict_hia(self, desc: Dict) -> float:
        """
        预测人肠道吸收
        基于分子量和TPSA
        """
        mw = desc.get('mw', 0)
        tpsa = desc.get('tpsa', 0)
        
        if mw > 500 or tpsa > 140:
            hia = 30 + (500 - mw) * 0.1 + (140 - tpsa) * 0.3
        else:
            hia = 90 - (mw - 200) * 0.02 - tpsa * 0.1
        
        return np.clip(hia, 0, 100)
    
    def _predict_pgp_substrate(self, desc: Dict) -> bool:
        """预测P-糖蛋白底物"""
        mw = desc.get('mw', 0)
        logp = desc.get('logp', 0)
        
        # 高分子量和高脂溶性更可能是P-gp底物
        return mw > 400 and logp > 3
    
    def _predict_vd(self, desc: Dict) -> float:
        """预测分布体积 (L/kg)"""
        logp = desc.get('logp', 0)
        
        # 高脂溶性药物分布体积较大
        vd = 0.5 + 0.3 * logp
        return np.clip(vd, 0.1, 5.0)
    
    def _predict_bbb(self, desc: Dict) -> float:
        """
        预测血脑屏障通透性
        LogBB = -0.0148*TPSA + 0.152*LogP + 0.139
        """
        tpsa = desc.get('tpsa', 0)
        logp = desc.get('logp', 0)
        
        logbb = -0.0148 * tpsa + 0.152 * logp + 0.139
        return logbb
    
    def _predict_ppb(self, desc: Dict) -> float:
        """预测血浆蛋白结合率 (%)"""
        logp = desc.get('logp', 0)
        
        # 高脂溶性药物蛋白结合率高
        ppb = 50 + 10 * logp
        return np.clip(ppb, 10, 99)
    
    def _predict_cyp_substrate(self, desc: Dict) -> Dict[str, bool]:
        """预测CYP酶底物"""
        mw = desc.get('mw', 0)
        logp = desc.get('logp', 0)
        
        results = {}
        for cyp in self.cyp_families:
            # 简化规则
            if cyp == '3A4':
                results[cyp] = mw > 300 and logp > 2
            elif cyp == '2D6':
                results[cyp] = desc.get('num_hacceptors', 0) > 2
            else:
                results[cyp] = False
        
        return results
    
    def _predict_cyp_inhibitor(self, desc: Dict) -> Dict[str, bool]:
        """预测CYP酶抑制剂"""
        logp = desc.get('logp', 0)
        
        results = {}
        for cyp in self.cyp_families:
            # 高脂溶性化合物更可能是CYP抑制剂
            results[cyp] = logp > 4
        
        return results
    
    def _predict_clearance(self, desc: Dict) -> float:
        """预测清除率 (mL/min/kg)"""
        mw = desc.get('mw', 0)
        logp = desc.get('logp', 0)
        
        # 低分子量和适中脂溶性清除快
        cl = 20 - (mw - 300) * 0.01 + (3 - abs(logp - 2)) * 2
        return np.clip(cl, 1, 50)
    
    def _predict_half_life(self, desc: Dict) -> float:
        """预测半衰期 (小时)"""
        mw = desc.get('mw', 0)
        
        # 高分子量药物半衰期长
        t_half = 2 + (mw - 300) * 0.01
        return np.clip(t_half, 0.5, 24)
    
    def _predict_ames(self, desc: Dict) -> bool:
        """预测Ames致突变性"""
        aromatic = desc.get('aromatic_rings', 0)
        
        # 多环芳烃可能致突变
        return aromatic > 3
    
    def _predict_herg(self, desc: Dict) -> bool:
        """预测hERG钾通道抑制"""
        logp = desc.get('logp', 0)
        tpsa = desc.get('tpsa', 0)
        
        # 高脂溶性和低TPSA易抑制hERG
        return logp > 4 and tpsa < 70
    
    def _predict_hepatotoxicity(self, desc: Dict) -> bool:
        """预测肝毒性"""
        mw = desc.get('mw', 0)
        logp = desc.get('logp', 0)
        
        # 大分子量和高脂溶性可能有肝毒性
        return mw > 500 and logp > 5
    
    def _predict_skin_sensitization(self, desc: Dict) -> bool:
        """预测皮肤敏感"""
        return False  # 简化模型
    
    def _predict_ld50(self, desc: Dict) -> float:
        """预测LD50 (mg/kg)"""
        mw = desc.get('mw', 0)
        logp = desc.get('logp', 0)
        
        # 经验公式
        ld50 = 1000 - 100 * abs(logp - 2) - (mw - 300) * 0.5
        return np.clip(ld50, 10, 5000)
    
    def _calculate_overall_score(self, caco2, hia, bbb, ames, herg, hepato, ld50) -> float:
        """计算ADMET综合评分"""
        score = 50.0
        
        # Absorption (+20分)
        if caco2 > -5.5:
            score += 10
        if hia > 70:
            score += 10
        
        # Distribution (+10分)
        if bbb > -1:
            score += 10
        
        # Toxicity (-30分)
        if ames:
            score -= 10
        if herg:
            score -= 10
        if hepato:
            score -= 10
        
        # LD50 (+10分)
        if ld50 > 500:
            score += 10
        
        return np.clip(score, 0, 100)
    
    def _get_default_profile(self) -> ADMETProfile:
        """返回默认ADMET Profile"""
        return ADMETProfile(
            caco2_permeability=0.0,
            hia=0.0,
            pgp_substrate=False,
            vd=0.0,
            bbb_permeability=0.0,
            ppb=0.0,
            cyp_substrate={cyp: False for cyp in self.cyp_families},
            cyp_inhibitor={cyp: False for cyp in self.cyp_families},
            clearance=0.0,
            half_life=0.0,
            ames_toxicity=True,
            herg_inhibition=True,
            hepatotoxicity=True,
            skin_sensitization=True,
            ld50=0.0,
            overall_score=0.0
        )


class DrugLikenessFilter:
    """类药性筛选器"""
    
    # PAINS子结构 (简化版)
    PAINS_SMARTS = [
        'C1(=O)C=CC(=O)C=C1',  # Quinone
        'c1ccc2c(c1)ncc(=O)[nH]2',  # Quinoline
        '[SH]',  # Thiol
        'N=N',  # Azo
        '[N+](=O)[O-]',  # Nitro
    ]
    
    def __init__(self):
        self.pains_patterns = [Chem.MolFromSmarts(s) for s in self.PAINS_SMARTS if Chem.MolFromSmarts(s)]
    
    def evaluate(self, mol: Chem.Mol) -> DrugLikenessProfile:
        """完整类药性评估"""
        try:
            # Lipinski
            lipinski_pass, lipinski_violations = self._check_lipinski(mol)
            
            # Veber
            veber_pass = self._check_veber(mol)
            
            # Ghose
            ghose_pass = self._check_ghose(mol)
            
            # Egan
            egan_pass = self._check_egan(mol)
            
            # Lead-like
            lead_like = self._check_lead_like(mol)
            
            # Fragment-like
            fragment_like = self._check_fragment_like(mol)
            
            # PAINS
            pains_alerts = self._check_pains(mol)
            
            # 综合评分
            overall = self._calculate_overall_score(
                lipinski_violations, veber_pass, pains_alerts
            )
            
            return DrugLikenessProfile(
                lipinski_pass=lipinski_pass,
                lipinski_violations=lipinski_violations,
                veber_pass=veber_pass,
                ghose_pass=ghose_pass,
                egan_pass=egan_pass,
                lead_like=lead_like,
                fragment_like=fragment_like,
                pains_alerts=pains_alerts,
                overall_score=overall
            )
            
        except Exception as e:
            logger.error(f"类药性评估失败: {e}")
            return self._get_default_profile()
    
    def _check_lipinski(self, mol: Chem.Mol) -> Tuple[bool, int]:
        """
        Lipinski五规则 (Ro5)
        MW <= 500, LogP <= 5, HBD <= 5, HBA <= 10
        """
        violations = 0
        
        try:
            mw = Descriptors.MolWt(mol)
            if mw > 500:
                violations += 1
        except:
            violations += 1
        
        try:
            logp = Crippen.MolLogP(mol)
            if logp > 5:
                violations += 1
        except:
            violations += 1
        
        try:
            hbd = Lipinski.NumHDonors(mol)
            if hbd > 5:
                violations += 1
        except:
            violations += 1
        
        try:
            hba = Lipinski.NumHAcceptors(mol)
            if hba > 10:
                violations += 1
        except:
            violations += 1
        
        return violations <= 1, violations
    
    def _check_veber(self, mol: Chem.Mol) -> bool:
        """
        Veber规则
        RotBonds <= 10, TPSA <= 140
        """
        try:
            rot_bonds = Lipinski.NumRotatableBonds(mol)
            tpsa = Descriptors.TPSA(mol)
            return rot_bonds <= 10 and tpsa <= 140
        except:
            return False
    
    def _check_ghose(self, mol: Chem.Mol) -> bool:
        """
        Ghose规则
        MW: 160-480, LogP: -0.4-5.6, Atoms: 20-70
        """
        try:
            mw = Descriptors.MolWt(mol)
            logp = Crippen.MolLogP(mol)
            atoms = mol.GetNumHeavyAtoms()
            
            return (160 <= mw <= 480 and 
                   -0.4 <= logp <= 5.6 and 
                   20 <= atoms <= 70)
        except:
            return False
    
    def _check_egan(self, mol: Chem.Mol) -> bool:
        """
        Egan规则
        LogP <= 5.88, TPSA <= 131.6
        """
        try:
            logp = Crippen.MolLogP(mol)
            tpsa = Descriptors.TPSA(mol)
            return logp <= 5.88 and tpsa <= 131.6
        except:
            return False
    
    def _check_lead_like(self, mol: Chem.Mol) -> bool:
        """
        先导化合物相似性
        MW: 250-350, LogP: 1-3
        """
        try:
            mw = Descriptors.MolWt(mol)
            logp = Crippen.MolLogP(mol)
            return 250 <= mw <= 350 and 1 <= logp <= 3
        except:
            return False
    
    def _check_fragment_like(self, mol: Chem.Mol) -> bool:
        """
        片段相似性
        MW: 120-250, LogP: -2-2.5
        """
        try:
            mw = Descriptors.MolWt(mol)
            logp = Crippen.MolLogP(mol)
            return 120 <= mw <= 250 and -2 <= logp <= 2.5
        except:
            return False
    
    def _check_pains(self, mol: Chem.Mol) -> int:
        """检查PAINS (泛活性化合物)"""
        alerts = 0
        try:
            for pattern in self.pains_patterns:
                if mol.HasSubstructMatch(pattern):
                    alerts += 1
        except:
            pass
        return alerts
    
    def _calculate_overall_score(self, lipinski_viol, veber_pass, pains) -> float:
        """计算综合评分"""
        score = 50.0
        
        # Lipinski (+30分)
        score += (4 - lipinski_viol) * 7.5
        
        # Veber (+20分)
        if veber_pass:
            score += 20
        
        # PAINS (-5分/警报)
        score -= pains * 5
        
        return np.clip(score, 0, 100)
    
    def _get_default_profile(self) -> DrugLikenessProfile:
        """返回默认Profile"""
        return DrugLikenessProfile(
            lipinski_pass=False,
            lipinski_violations=4,
            veber_pass=False,
            ghose_pass=False,
            egan_pass=False,
            lead_like=False,
            fragment_like=False,
            pains_alerts=99,
            overall_score=0.0
        )


if __name__ == '__main__':
    # 测试
    from rdkit import Chem
    
    # 测试分子: Aspirin
    mol = Chem.MolFromSmiles('CC(=O)Oc1ccccc1C(=O)O')
    
    # ADMET预测
    admet_predictor = ADMETPredictor()
    admet = admet_predictor.predict(mol)
    
    print("ADMET评估:")
    print(f"  Caco-2: {admet.caco2_permeability:.2f}")
    print(f"  HIA: {admet.hia:.1f}%")
    print(f"  BBB: {admet.bbb_permeability:.2f}")
    print(f"  Ames: {admet.ames_toxicity}")
    print(f"  hERG: {admet.herg_inhibition}")
    print(f"  综合评分: {admet.overall_score:.1f}/100")
    
    # 类药性评估
    druglikeness = DrugLikenessFilter()
    profile = druglikeness.evaluate(mol)
    
    print("\n类药性评估:")
    print(f"  Lipinski: {profile.lipinski_pass} (违规: {profile.lipinski_violations})")
    print(f"  Veber: {profile.veber_pass}")
    print(f"  PAINS: {profile.pains_alerts} 警报")
    print(f"  综合评分: {profile.overall_score:.1f}/100")
