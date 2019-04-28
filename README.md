A. Resources:
- README.md: this file.
- GSE18123GPL570.csv and GSE18123GPL6244.csv: two datasets
- result.csv: result for each method with different number of features. Measured by accuracy, auc, f1, and	mcc.

- source codes:
+ conventional_feature_ranking.py: 12 conventional feature ranking methods plus AUC (baseline) and SFARI.
+ gene_interaction_co_expression.py: create co-expression networks.
+ gene_interaction_pathway.py: create pathway networks ('DO', 'GObp','GOcc', 'GOmf', and 'HPO' networks).
+ gene_interaction_PPI.py: create PPI networks.
+ gene_interaction_pubmed.py: create pubmed-embedding networks.
+ gene_embedding.py: create word, hence gene name included, embedding, used as input for gene_interaction_pubmed.py
+ gene_interaction_joint.py: create a ensembling network where an edge is defined between two genes if it is existed in more than half of all the 8 networks ('co_expression','DO', 'GObp','GOcc', 'GOmf', 'HPO','PPI', and 'pubmed').
+ geneRank.m: score genes based on their networks: co-expression, pathway, PPI, and pubmed-embedding networks
+ score_2_index.py: turn the score returned by geneRank.m into the index of genes in the datasets. 
+ python classify.py: running the classification task with the features ranked highest by the above methods, be it conventional or network-based.

- folder pubmed/: 
+ pubmed.txt (unzip pubmed.7z): contain PubMed raw text as input to learn word/gene embeddings.
+ gene_embedding.embed (unzip gene_embedding.7z): word/gene embeddings returned by running word2vec models on pubmed.txt.
+ gene_name_ref.txt: mapping vocabulary (gene names) between gene_embedding.embed and the datasets. 

- folder mapping/ (unzip mapping.7z):
+ hgncTo*.csv: ontology mapping from genes to DO, GO, and HPO.
+ hippie_current.txt: protein-protein interactions, downloaded from http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/hippie_current.txt 
+ SFARI-Gene_genes_export09-06-2018.csv: SFARI list of genes, downloaded from https://gene.sfari.org//wp-content/themes/sfari-gene/utilities/download-csv.php?api-endpoint=genes

- gene_interaction/ (upzip all *.7z files into *.csv files): store the output of gene_interaction_*.py, creating gene-gene networks. 
+ gene_interaction/score/: store the output of geneRank.m, scoring genes by their networks.

- feat_ranking/ (unzip feat_ranking.7z): gene ranking, either return by conventional_feature_ranking.py or score_2_index.py, which is the network-based feature ranking.
This is also the input for classify.py, which need a list of features for their classification.    

B. Step-by-step running:

0. Installing Python libaries needed for 
- conventional feature selections:
pip install skfeature-chappers
- training embbeded vectors: 
pip install gensim

1. Running conventional feature ranking:
python conventional_feature_ranking.py index

where index is in the range of the list of conventional methods ['chi_square','cmim', 'f_score','fisher_score','gini_index','icap','jmi','ll_l21','ls_l21','reliefF','rfs','trace_ratio', 'SFARI','AUC']
The material for SFARI is stored at /home/thinng/code/2019/bib/mapping/SFARI-Gene_genes_export09-06-2018.csv, 
downloaded from https://gene.sfari.org//wp-content/themes/sfari-gene/utilities/download-csv.php?api-endpoint=genes
The ranking of genes will be stored at the folder of feat_ranking/

2. Building gene-gene networks:

- co_expression interactions:
python gene_interaction_co_expression.py

- pathway interactions:
python gene_interaction_pathway.py

The mapping of gene-ontologies is stored at mapping/hgncTo*.csv. There are five ontologies: 'DO', 'GObp','GOcc', 'GOmf', and 'HPO'.

- PPI networks:
python gene_interaction_PPI.py

The material for building PPI is stored at mapping/hippie_current.txt, downloaded from http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/hippie_current.txt

- PubMed based networks:
For learning gene-gene interactions from PubMed, first the embedding for genes should be learned by running:
python gene_embedding.py

Raw data extracted from PubMed is stored at /pubmed/pubmed.txt. The embbed model is stored at pubmed/gene_embedding.embed.
Then the network is learned through:
python gene_interaction_pubmed.py

- Joint networks:
For ensembling network-based features. An edge is defined between two genes if it is existed in more than half of all the networks.
There are 8 networks in total: 'co_expression','DO', 'GObp','GOcc', 'GOmf', 'HPO','PPI', and 'pubmed'. 
python gene_interaction_joint.py

- Running all these python programs result in gene-gene networks stored as gene_interaction/

3. Running network-based feature ranking for these networks by:
geneRank.m

The score for genes returned by running geneRank.m is stored at /gene_interaction/score/. We convert the score to index for genes by running
python score_2_index.py

The ranking of genes returned by running score_2_index.py will be stored at the folder of feat_ranking/

4. Now we have feature ranking for both conventional and network-based methods. Then the features with highest ranks will be used in the following classification. 
python classify.py

The result for all these feature selections returned by running classify.py will be stored in result.csv.

