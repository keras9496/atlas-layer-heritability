"""Lineage assignment rules of experiment 3c (PROTOCOL_exp3c_2026-09-28), extracted verbatim from z01_lineages.py.
 The original script also retrieved track metadata through the AlphaGenome API (scorer_metadata); that part is not included.
 assign(r) takes a metadata row with biosample_type and biosample_name and returns the lineage, "excluded_cell_line" or "unassigned".
 Rules are applied in order; the first match wins (case-insensitive). See Supplementary Table S2 for the resulting track counts.
"""
import re

EXCLUDE = r"^blood$|venous blood|peripheral blood mononuclear|^mononuclear cell|bone marrow|umbilical cord blood|cerebrospinal fluid|adult organism|forelimb$|hindlimb$|mesendoderm|mesodermal cell|endodermal cell|ectodermal cell|cell line|hematopoietic multipotent|reticulocyte"
RULES = [  # (계통, 정규식) — 순서가 중요하다
    ("CNS_tissue", r"TISSUE::(?!.*(renal|kidney)).*(brain|cortex|gyrus|lobe$|occipital pole|caudate|putamen|globus pallidus|nucleus accumbens|amygdala|hippocamp|ammon|insula|thalamus|diencephalon|cerebell|pons|medulla oblongata|substantia nigra|locus ceruleus|corpus callosum|spinal cord|germinal matrix)"),
    ("neural_cells", r"CELL::.*(\bneuron|\bneural|neurosphere|astrocyte|oligodendrocyte|schwann|purkinje)"),
    ("myeloid", r"CELL::.*(monocyte|macrophage|dendritic|langerhans|neutrophil|eosinophil|basophil|mast cell|myelocyte|common myeloid progenitor)"),
    ("lymphoid", r"CELL::.*(\bt cell|\bt-cell|\bb cell|natural killer|t-helper|t follicular|lymphoblast)|TISSUE::.*(thymus|spleen|tonsil|lymph node|peyer)"),
    ("liver", r"hepatocyte|TISSUE::.*liver"),
    ("GI", r"TISSUE::.*(esophag|stomach|gastroesophageal|colon|colonic|intestin|duoden|rectum|rectal|appendix|pancrea|gallbladder|salivary|saliva|parotid|submandibular|tongue|mouth|lingual gland|fungiform)|CELL::.*(intestinal epithelial|epithelial cell of esophagus|acinar cell of salivary|islet|pancrea|mouth mucosa|smooth muscle cell of colon|smooth muscle cell of the esophagus|enteric smooth muscle)"),
    ("cardiovascular", r"TISSUE::.*(heart|myocardium|atri|ventricle|cardiac septum|valve|aorta|artery|vein|vena cava)|CELL::.*(cardiac muscle|cardiac myocyte|endothelial|pericyte|aortic smooth|smooth muscle cell of the .*(arter|vascul|vein)|mesothelial cell of epicardium|epithelial cell of umbilical artery)"),
    ("skeletal_muscle", r"TISSUE::(?!smooth).*(muscle|gastrocnemius|psoas|soleus|diaphragm)|CELL::.*(skeletal muscle|myoblast|myotube|myocyte|satellite cell)"),
    ("adipose", r"adipo|fat cell|fat pad|preadipocyte"),
    ("stromal", r"CELL::.*(fibroblast|mesenchymal|chondrocyte|osteoblast|osteocyte|tendon|synovial|stromal|stellate|smooth muscle|myometrial|dermal papilla|keratocyte|annulus pulposus|nucleus pulposus|trabecular meshwork|leptomeningeal|perineural)|TISSUE::.*(tendon|ligament|femur|dura mater|meninx|smooth muscle tissue)"),
    ("reproductive_endocrine", r"TISSUE::.*(testis|ovary|gonad|\buter|cervix|vagina|placenta|chorion|amnion|trophoblast|epididymis|seminal|vas deferens|penis|adrenal|thyroid|parathyroid|pituitary|pineal|umbilical cord|fallopian)|CELL::.*(trophoblast|sertoli|ovary|amnio|chorion|placental)"),
    ("epithelial_barrier", r"keratinocyte|melanocyte|epiderm|skin|outer root sheath|CELL::.*(bronchial epithelial|tracheal epithelial|respiratory epithelial|airway epithelial|alveol|nasal|pneumocyte|\bkidney|\brenal|glomerular|nephron|proximal tubule|urothelial|mammary|prostate|epithelial cell of malassez|gingival|corneal|lens|retinal pigment|iris|ciliary|olfactory|choroid plexus|mesangial|sebaceous)|TISSUE::.*(lung|trachea|\bkidney|\brenal|metanephros|ureter|urethra|bladder|breast|prostate|eye|retina|vitreous|nail|throat|olfactory|\bmole\b)"),
]

def assign(r):
    if r.biosample_type == "cell_line": return "excluded_cell_line"
    s = f"{'TISSUE' if r.biosample_type == 'tissue' else 'CELL'}::{r.biosample_name}"
    if re.search(EXCLUDE, r.biosample_name, re.I): return "unassigned"
    for g, p in RULES:
        if re.search(p, s, re.I): return g
    return "unassigned"
