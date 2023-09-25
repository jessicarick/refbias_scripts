# simulation schematic
# figure 1 in rick et al.

#devtools::install_github('rich-iannone/DiagrammeR')
install.packages("DiagrammeRsvg")
install.packages("rsvg")

library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)


DiagrammeR::grViz("
digraph boxes_and_circles {

  # a 'graph' statement
  graph [overlap = true, fontsize = 14]

  # several 'node' statements
  node [shape=plaintext,
        fontname = Helvetica,
        penwidth=3,
        fontcolor='#0A9396',
        fontsize=18]
  s1 [label = 'SimPhy'];
  s2 [label = 'TreeToReads'];
  s3 [label = 'bwa'];
  s4 [label = 'samtools\nbcftools\nvcftools'];
  s5 [label = 'RAxML\nASTRAL-III'];
  s6 [label = 'R'];
  
  node [shape = box,
        fontname = Helvetica,
        #fixedsize=true,
        width=4,
        color='#0A9396',
        penwidth=3,
        style=filled,
        fillcolor='#0A9396',
        fontcolor=white]
  a [label = 'Simulate 10 species trees\n(50 tip taxa) for\neach tree height',fontsize=18]; 
  
  node [shape = box,
        fontname = Helvetica,
        #fixedsize=true,
        width=3,
        nodesep=1,
        color='#0A9396',
        penwidth=3,
        fillcolor=white,
        fontcolor=black,
        fontsize=16]
  b [label = 'Simulate gene geneologies\n for 2000 loci constrained to each\nspecies tree',width=4];
  bc [label = 'Evolve gene sequences along gene genealogies;\ndraw variable sites from normal distribution N(150,30)',width=5.5];
  c [label = 'Size-select ~500bp fragments and simulate\n150bp fastq reads for each locus for each individual',width=4];
  
  e [label = 'Align fastq reads to\neach reference genome'];
  f [label = 'Call variants and filter VCF for\nQUAL=40, biallelic SNPs,\nMAC, and missing data', width=4];
  h [label = 'Remove invariant sites,\nconvert to phylip'];
  i [label = 'Maximum likelihood\ntree inference'];
  i2 [label = 'Split by locus; infer\ngene trees'];
  i3 [label = 'Species tree inference'];
  j [label = 'Summary statistics and modeling\nfor all trees for each reference genome\nand analysis method choice\nfrom each of 10 simulations']

  node [shape = oval,
        fontname = Helvetica,
        nodesep=0.25,
        style=dashed,
        color=gray]
  a1 [label = 'Speciation rates:\n0.00001 (high ILS)\n0.000005 (med ILS)\n0.000001 (low ILS)'];
  d [label = 'Reference genomes:\noutgroup taxon (EXT),\nrandom ingroup taxon (INT)'];

  g3 [label = 'Min minor allele count (MAC):\n0 / 1 / 2 / 3 / 4 / 5 / 10\n\nMissing data allowed:\n0 / 0.25 / 0.50 / 0.75 / 0.90'];





  # 'edge' statements
  edge [color = slategrey]
  a->b
  b->bc->c
  c->e->f->h
  h->i
  h->i2->i3->j
  i->j

  edge [color = slategrey,arrowhead=none]
  #a->a1 [len=0.5]
  #a1->b [minlen=1]
  #g4->f
  #g3->f
  #g3->g4
  #e->d->f
                  { rank = same; d;e }
                  { rank = same; f;g3 }
                  { rank = same; s1;a }
                  { rank = same; s2;bc }
                  { rank = same; s3;e }
                  { rank = same; s4;f }
                  { rank = same; s5;i;i3 }
                  { rank = same; s6;j }

}
") %>%
  export_svg %>%
  charToRaw %>%
  rsvg_pdf("~/Downloads/simulation_schematic.pdf")

graph %>% DiagrammeR::export_graph(file_name = "~/Downloads/simulation_schematic.pdf")
