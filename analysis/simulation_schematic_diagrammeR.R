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
  a [label = 'Simulate species tree\n(50 tip taxa)',fontsize=18]; 
  
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
  b [label = 'Simulate gene geneologies\n for 2000 loci (5000bp each)\nconstrained to each species tree',width=4];
  bc [label = 'Draw variable sites\nfrom N(150,30)',width=3];
  c [label = 'Simulate 150bp fastq reads for\neach individual at each locus',width=4];
  
  e [label = 'Align fastq reads to\nreference genomes'];
  f [label = 'Call variants and filter VCF for\nQUAL=40, biallelic SNPs only', width=4];
  h [label = 'Remove invariant sites,\nconvert to phylip'];
  i [label = 'Maximum likelihood\ntree estimation (RAxML)'];
  j [label = 'Summary statistics\nand analysis']

  node [shape = oval,
        fontname = Helvetica,
        nodesep=0.25,
        style=dashed,
        color=gray]
  a1 [label = 'Speciation rate\n0.00001 (high ILS)\n0.000005 (med ILS)\n0.000001 (low ILS)'];
  d [label = 'Choose outgroup taxon and\none random ingroup individual\nas reference genomes'];

  g3 [label = 'Minor allele count\n0 / 1 / 2 / 3 / 4 / 5 / 10'];
  g4 [label = 'Missing data\n0 / 0.25 / 0.50 / 0.75 / 0.90'];

  # several 'edge' statements
  edge [color = slategrey]
  b->bc->c
  c->d->e->f->g3->g4->h->i
  i->j

  edge [color = slategrey,arrowhead=none]
  a->a1 [len=0.5]
  a1->b[minlen=1]
}
")
