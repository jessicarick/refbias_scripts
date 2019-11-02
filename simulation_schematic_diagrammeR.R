devtools::install_github('rich-iannone/DiagrammeR')

grViz("
digraph boxes_and_circles {

  # a 'graph' statement
  graph [overlap = true, fontsize = 14]

  # several 'node' statements
  node [shape = box,
        fontname = Helvetica,
        #fixedsize=true,
        width=3,
        color=midnightblue,
        penwidth=3,
        style=filled,
        fillcolor=midnightblue,
        fontcolor=white]
  a [label = 'Simulate species tree\n(100 individuals)',fontsize=16]; 
  
  node [shape = box,
        fontname = Helvetica,
        #fixedsize=true,
        width=2.5,
        nodesep=1,
        color=midnightblue,
        penwidth=3,
        fillcolor=white,
        fontcolor=black]
  b [label = 'Simulate 1000 gene\ngeneologies for 1000bp loci\nconstrained to species tree',width=3];
  c [label = 'Simulate 150bp fastq reads\nfor each individual at\neach locus with\nrandom mutation rate',width=3]; 
  
  e [label = 'Align fastq reads to\nreference genome'];
  e2 [label = 'Align fastq reads to\nreference genome'];
  f [label = 'Call variants'];
  f2 [label = 'Call variants'];
  g [label = 'Filter VCF\nQUAL=40,\nBiallelic SNPs only'];
  g2 [label = 'Filter VCF\nQUAL=40,\nBialleleic SNPs only'];
  h [label = 'Remove invariant sites,\nconvert to phylip'];
  h2 [label = 'Remove invariant sites,\nconvert to phylip'];
  i [label = 'Maximum likelihood\ntree estimation (RAxML)'];
  j [label = 'Summary statistics\nand analysis']

  node [shape = oval,
        fontname = Helvetica,
        nodesep=0.25,
        style=dashed,
        color=gray]
  a1 [label = 'Tree height\n500,000\n2,000,000\n10,000,000'];
  d [label = 'Choose random ingroup\nindividual as reference'];
  d2 [label = 'Choose outgroup\nindividual as reference'];

  g3 [label = 'Minor allele frequency\n0 / 0.01 / 0.02 / 0.03 / 0.04 / 0.05 / 0.10'];
  g4 [label = 'Missing data\n0 / 0.25 / 0.50 / 0.75 / 1.0'];

  # several 'edge' statements
  edge [color = slategrey]
  b->c
  c->d->e->f->g->{g3,g4}->h->i
  c->d2->e2->f2->g2->{g3,g4}->h2->i
  i->j

  edge [color = slategrey,arrowhead=none]
  a->a1 [len=0.5]
  a1->b[minlen=2]
}
")
