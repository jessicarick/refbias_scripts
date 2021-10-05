#devtools::install_github('rich-iannone/DiagrammeR')

DiagrammeR::grViz("
digraph boxes_and_circles {

  # a 'graph' statement
  graph [overlap = true, fontsize = 14]

  # several 'node' statements
  node [shape = box,
        fontname = Helvetica,
        #fixedsize=true,
        width=4,
        color='#0A9396',
        penwidth=3,
        style=filled,
        fillcolor='#0A9396',
        fontcolor=white]
  a [label = 'Simulate species tree\n(50 tips)',fontsize=18]; 
  
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
  b [label = 'Simulate 2000 gene\ngeneologies for 5000bp loci\nconstrained to species tree',width=3];
  c [label = 'Simulate 150bp fastq reads\nfor each individual at\neach locus with\nrandom mutation rate',width=3]; 
  
  e [label = 'Align fastq reads to\nreference genome'];
  f [label = 'Call variants\nand filter VCF for\nQUAL=40,\nbiallelic SNPs only'];
  h [label = 'Remove invariant sites,\nconvert to phylip'];
  i [label = 'Maximum likelihood\ntree estimation (RAxML)'];
  j [label = 'Summary statistics\nand analysis']

  node [shape = oval,
        fontname = Helvetica,
        nodesep=0.25,
        style=dashed,
        color=gray]
  a1 [label = 'Speciation rate\n0.00001 (high ILS)\n0.000005 (med ILS)\n0.000001 (low ILS)'];
  d [label = 'Choose random ingroup\nindividual as reference'];
  d2 [label = 'Choose outgroup\nindividual as reference'];

  g3 [label = 'Minor allele count\n0 / 1 / 2 / 3 / 4 / 5 / 10'];
  g4 [label = 'Missing data\n0 / 0.25 / 0.50 / 0.75 / 0.90'];

  # several 'edge' statements
  edge [color = slategrey]
  b->c
  c->{d,d2}->e->f->g3->g4->h->i
  i->j

  edge [color = slategrey,arrowhead=none]
  a->a1 [len=0.5]
  a1->b[minlen=1]
}
")
