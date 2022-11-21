# REFMAKER: create your own reference to target nuclear sequences in shotgun libraries

REFMAKER is a user-friendly pipeline providing different tools to create nuclear references genomic assemblies of shotgun libraries.

REFMAKER is a command-line program that needs to be run from a terminal/console, by calling different tasks, called 'modes', along with another parameter corresponding to specific 'options' (see Figure 1). It can be parameterized in order to:
 1. perform the assembly for each library
 2. perform the meta-assembly of these libraries to create a catalog
 3. clean the catalog by selecting the nuclear regions, and by removing duplicates from clustering steps
 4. map the libraries into the cleaned catalog
 5. call the variants
 6. get the consensus sequences
 7. filter the final sequences to create a phylogenetic matrix



 <b>REFMAKER flowchart</b>


 ![Fig.1. REFMAKER worflow](./resources/img/Refmaker_workflow.png)
 >**Fig. 1. REFMAKER workflow**.
