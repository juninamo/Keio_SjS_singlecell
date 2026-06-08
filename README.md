![](https://komarev.com/ghpvc/?username=juninamo&style=flat-square&color=green&label=REPOSITORY+VIEWS)


# Multi-modal single-cell analysis of salivary glands from Sjögren's disease by different autoantibody profiles
Analytic codes for multi-omics data (scRNA-seq, TCR/BCR repertoire, and spatial transcriptome) from patients with Sjögren's disease, Consortium for Drug Discovery for Immuno-Inflammatory Diseases

[![Article DOI](https://img.shields.io/badge/Nature%20Communications-10.1038%2Fs41467--025--63935--9-blue)](https://doi.org/10.1038/s41467-025-63935-9)
[![Data DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20588865.svg)](https://doi.org/10.5281/zenodo.20588865)

&nbsp;&nbsp;


## Comparative single-cell and spatial profiling of anti-SSA-positive and anti-centromere-positive Sjögren's disease reveals common and distinct immune activation and fibroblast-mediated inflammation

Sjögren's disease (SjD) is an autoimmune disease that causes salivary gland dysfunction due to immune-mediated destruction. While autoantibodies such as anti-SSA and anti-centromere (CENT) are associated with distinct clinical manifestations, the molecular features remain to be elucidated. In this study, we apply multi-modal single-cell technologies: single-cell RNA sequencing, T cell and B cell receptor sequencing and spatial transcriptomics to salivary gland lesions, aiming to elucidate common and unique cellular and transcriptional signatures linked to different autoantibody profiles. Our analysis demonstrates that GZMB⁺GNLY⁺ CD8⁺ T cells are the main expanded subset across different autoantibody statuses, highlighting their central role in SjD pathogenesis, while the enrichment of memory B cells is more prominent in anti-CENT-positive patients. Cytokine signaling also differs by autoantibody profile, with an activated interferon signature in anti-SSA-positive patients, whereas TGFβ signaling is enhanced in anti-CENT-positive patients. Furthermore, spatial profiling reveals THY1⁺ fibroblasts, expressing complement genes and chemokines, as key hubs orchestrating inflammation within the salivary glands. These findings deepen our understanding of the pathogenesis of SjD, and may inform the development of targeted and personalized therapeutic strategies.

&nbsp;&nbsp;

## Citation

Inamo J\*, Takeshita M\*, Suzuki K, Tsunoda K, Usuda S, Kuramoto J, Moody J, Hon CC, Ando Y, Sasaki T, Yoshitake K, Mitsuyama S, Asakawa S, Kanai Y, Takeuchi T, Kaneko Y. **Comparative single-cell and spatial profiling of anti-SSA-positive and anti-centromere-positive Sjögren's disease reveals common and distinct immune activation and fibroblast-mediated inflammation.** *Nature Communications* **16**, 8299 (2025). doi:[10.1038/s41467-025-63935-9](https://doi.org/10.1038/s41467-025-63935-9). PMID: [40983612](https://pubmed.ncbi.nlm.nih.gov/40983612/).

\* These authors contributed equally.

&nbsp;&nbsp;

## Data availability

The processed data supporting this study are openly available on Zenodo:
**[10.5281/zenodo.20588865](https://doi.org/10.5281/zenodo.20588865)**.

The deposit includes:

- scRNA-seq raw count matrix (all genes, sparse 10x Matrix Market format) and per-cell metadata (sample, clinical and cell-type annotations)
- Single-cell TCR and BCR repertoire tables (AIRR rearrangement format)
- Spatial transcriptomics (10x Visium) per-slide count matrices and spot coordinates
- Clinical metadata for the Visium samples

&nbsp;&nbsp;

<!--
<kbd>
<img src="https://github.com/juninamo/consortium_step2_SjS/blob/main/images/Figure1_overview.png" width="800" align="center">
</kbd>

&nbsp;&nbsp;

**Figure 1. Study design**
-->

&nbsp;&nbsp;

## Contact us
Please contact us (Jun Inamo: juninamo@keio.jp) with any questions or comments.

The data presented here comes from the [Division of Rheumatology, Department of Internal Medicine, Keio University School of Medicine](https://www.med.keio.ac.jp/en/) through collaborating with [RIKEN](https://www.riken.jp/en/) and University of Tokyo.

<!--
<kbd>
<img src="https://github.com/juninamo/consortium_step2_SjS/blob/main/images/Keio_logo.png" width="300" align="center">
</kbd>
-->

## Acknowledgments
This research was performed in the Immune-mediated Inflammatory Diseases Consortium for Drug Development, supported by Keio University School of Medicine, Division of Rheumatology and Gastroenterology, Iwate Medical University,Division of Allery and Rheumatology, Department of Internal Medicine. the National Institute of Biomedical Innovation (NIBIO), Toxicogenomics Informatics Project, Ono Pharmaceutical Co., Ltd., Daiichi Sankyo Co., Ltd., Mitsubishi Tanabe Pharma CO.. This publication is part of the Single Cell Medical Network of Japan. This research was also funded by Japan Society for the Promotion of Science (JSPS) Grant-in-Aids for JSPS Fellows (grant number JP21J00596), JSPS KAKENHI (Grant Number JP20K17430, JP20H03720, and JP22K08528), Keio University Academic Development Funds (M.T.), The Mochida Memorial Foundation for Medical and Pharmaceutical Research, Kowa Life Science Foundation, Japan Rheumatism Foundation, Takeda Science Foundation. 

&nbsp;&nbsp;
