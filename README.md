#  GRIMM: Graph IMputation and Matching for HLA Genotypes

Graph-based imputation and matching, a.k.a. new imputation method (NIM).

| Subdirectory    | Description |
| :-------------- | :---------- |
| graph_generator | Python code to generate (haplotype frequency) CSV files which can be imported into Neo4j |
| gl_impute       | old Python code for graph-based imputation.  Still needed for cypher queries (?) |
| multi_race_ipute | Python code for graph-based imputation |
| validation      | Code for validation of components of this repository |

# Walkthrough

### 1. Generate Initial Graph
 Follow [Graph Generation](graph_generator/README.md) to generate haplotype frequency CSV files and import to Neo4j.


### 2. Perform Imputation
 Follow [Imputation](multi_race_impute/README.md) directions to produce imputation results from frequency files.


### 3. Perform Matching
  Follow [Match Validation](matching/README.md) to perform matching and compare against consensus results.
