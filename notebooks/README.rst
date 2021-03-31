Directory for the Jupyter notebooks.

Naming convention is a number (for ordering), the creator's initials, and a short ``-`` delimited description, e.g. ``1.0-zbh-initial-data-exploration``.

Notebooks
-----------------------


1.0-ma-uniprot-data-extraction.ipynb
  Initial attempt to extract data from Uniprot database (using random enzymes)
1.2-ma--sp--uniprot-data-extraction-optimized :
  From Excel file of UniProt ID's, extracts ID and forms url and forms soup
  extracts Km values and Vm values and inputs them into separate Panda dataframe
1.2.1-sp--ma--uniprot-data-extraction-optimized
  Attaches UniProt ID's as well and forms an excel table at the end
1.2.5-sp--ma--uniprot-data-extraction-optimized
  LATEST : Combined Km and Vm values, forms data table at the end with split values
1.4-ma-pathway-notebook-creation-workflow.ipynb
  Initial attempt to create a pathway construction workflow using Amino Acid reactions (does not work) 
1.0-sz--rm-Keq_Extraction_for_Reactions.ipynb
  From Excel file of reaction(rxn) BIGG IDs for all pathways, extract rxn IDs to form BIGG urls.
  Then, extracts EC numbers for each rxn from BIGG website.
  Use the EC numbers to extract all Keq datas from Sabio-RK website and map the EC numbers result with rxn IDs.
1.1-ma-Keq_eQuilibrator.ipynb
  Same details as 1.0-sz--rm-Keq_Extraction_for_Reactions.ipynb
  Updated to include Beautiful Soup code to extract exact Keq data 
2.0-sz--rm-enzyme-data-extraction-sabiork.ipynb
  (Final Version for Extraction)
  From Excel file of "EnzymeModulesList", extract UniProt IDs for enzymes to form Sabio-RK urls.
  Use Helium to extract enzyme data from Sabio-RK website
1.0-sz--rm--Keq_eQuilibrator.ipynb
  From Excel file of "rxn_testing", extract reaction BIGG IDs to form BIGG urls in order to get full reaction names and the first EC number.
  Use the EC numbers to form eQuilibrator urls. Automate to adjust pH, pMg, and ionic strength (pseudo parameter for testing; still need to find standards)
  May need to manually extract the Gibbs free energy (no success with Helium and BeautifulSoup yet)
 
  
