Important notes on this pipeline:
1. We are using the TOGA-based annotations for these genomes.  They are, in some cases (potentially due to massive many-to-many orthology) horrifically wrong.  
For example, the Tursiops_truncatus gtf claims that the CDS for ZNF577 (an ~800 amino acid protein) is nearly 200,000 nucleotides long.
It might be a good idea to remove genes with badly defined orthology somehow (or just remove genes with way more nucleotides assigned to them than their known length in humans).
It makes sense to remove these genes after running the test, this is why we also back up the final bed file of CDS (for comparison to human gene lengths later)