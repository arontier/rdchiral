from rdkit import Chem
from rdkit.Chem import AllChem
from rdchiral.main import rdchiralReaction, rdchiralReactants, rdchiralRun

from rdkit import Chem
from rdchiral.main import rdchiralRunText

from rdchiral.main import rdchiralRunText

# 매핑/label 없이 단순 반응 SMARTS로 먼저 테스트
rxn_smarts = "CCCBr.[OH-]>>CCCO.[Br-]"   # 아주 단순 패턴
reactant_smiles = "CCCBr.[OH-]"

out_plain, out_mapped = rdchiralRunText(
    rxn_smarts,
    reactant_smiles,
    return_mapped=True
)

print("plain:", out_plain)
print("mapped:", out_mapped)


rxn_smarts = "[CH3:1][CH2:2][Br:3].[OH-:4]>>[CH3:1][CH2:2][O:4].[Br-:3]"

mol1 = Chem.MolFromSmiles("CCCBr")
mol2 = Chem.MolFromSmiles("[OH-]")

smis = [Chem.MolToSmiles(m, isomericSmiles=True) for m in (mol1, mol2)]
print("Individual SMILES:", smis)

reactant_smiles = ".".join(smis)
print("reactant_smiles:", reactant_smiles)
print("Num fragments in reactant_smiles:", len(reactant_smiles.split(".")))

# 여기서부터는 RDChiral이 알아서 rdchiralReaction/rdchiralReactants를 내부에서 생성
out_plain, out_mapped = rdchiralRunText(
    rxn_smarts,
    reactant_smiles,
    return_mapped=True
)

print("=== RDChiral 결과 ===")
for p, m in zip(out_plain, out_mapped):
    print("plain :", p)
    print("mapped:", m)
    print("-----")

# 1. RDKit 쪽: 이미 Mol을 가지고 있다고 가정
# (예시로 여기서 SMILES→Mol을 만들지만, 보통은 앞 단계에서 만든 mol을 사용할 것)
reactant_mols = [
    Chem.MolFromSmiles("CCBr"),   # 예: 1-bromopropane
    Chem.MolFromSmiles("[OH-]")   # hydroxide
]

# RDChiral은 atom-mapped를 권장하지만, 여기서는 단순 예시로 비-mapped 사용
# 실제 템플릿 적용/역합성에서는 atom mapping을 반드시 넣어야 함[web:64]
reactant_smiles = ".".join(Chem.MolToSmiles(mol, isomericSmiles=True)
                           for mol in reactant_mols)
print(reactant_smiles)

# 2. RDChiral 쪽: reaction SMARTS 정의
# Correct SN2: [OH-] + alkyl-Br -> alcohol (neutral -O-H after protonation) + Br-
rxn_smarts = "[CH3:1][CH2:2][Br:3].[OH-:4]>>[CH3:1][CH2:2][O:4].[Br-:3]"

# RDChiral 객체들 초기화
rxn = rdchiralReaction(rxn_smarts)                # RDKit ChemicalReaction 래핑
reactants = rdchiralReactants(reactant_smiles)    # RDKit Mol + achiral 복사본 생성[web:60][web:67]

# 3. RDChiral로 반응 실행 (Mol 기반 파이프라인에 통합)
outcomes_plain, outcomes_mapped = rdchiralRun(
    rxn,
    reactants,
    return_mapped=True     # mapped 결과도 같이 받기
)

print("=== RDChiral 결과 ===")
for plain, mapped in zip(outcomes_plain, outcomes_mapped):
    print("plain :", plain)
    print("mapped:", mapped)
    print("-----")

# 4. RDKit Mol로 다시 가져오고 싶다면:
product_mols = [Chem.MolFromSmiles(s) for s in outcomes_plain]