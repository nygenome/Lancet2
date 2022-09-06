# Visualization

By passing the ```--graphs-dir``` parameter and providing a directory where you want the graphs to be placed, Lancet will write deBruijn graphs for each window inspected and place them into the directory as dot files. NOTE: whichever directory is given as the graphs-dir will be cleared so be mindful of what directory you provide.

```bash
./lancet2 pipeline -t tumor.bam -n normal.bam -r ref.fasta --graphs-dir ./dot_graphs_dir
```

The above command will export the DeBruijn graph at various stages of the assembly (low coverage removal, graph compression and traversal) to the following set of files:

1. chr:start-end_c0_raw_graph.dot
2. chr:start-end_cX_before_compression.dot
3. chr:start-end_cX_after_compression.dot
4. chr:start-end_cX_path_flow.dot

Where X refers to the connected component within the graph (in most cases only one).

These files can be rendered using the dot utility available in the [Graphviz](http://www.graphviz.org/) visualization software package.

```
dot -Tpdf -o example_file.pdf example_file.dot
```

The above command will create a example_file.pdf file that shows the graph. For large graphs, Adobe Acrobat Reader may have troubles rendering the graph in which case we recommend opening the PDF file using the "Preview" image viewer software available in MacOS.

Below is an example of what the generated graphs may look like. The blue nodes are k-mers shared by both tumor and normal; the white nodes are k-mer with low support (likely sequencing errors); the green nodes are k-mers only present in the normal; the red nodes are k-mers only present in the tumor.:

The first image below shows the raw graph before removing the low coverage nodes.
![raw_graph](https://github.dev/nygenome/Lancet2/blob/db225350de4da2a125694e03e01dbc006a9865fc/website/static/img/chr14_72547800-72548098_c0_raw_graph.png)

The second image below shows the graph before compression and tip removal.
![before_compression](https://github.dev/nygenome/Lancet2/blob/db225350de4da2a125694e03e01dbc006a9865fc/website/static/img/chr14_72547800-72548098_c1_before_compression.png)

The third image below shows the graph after compression and tip removal.
![after_compression](https://github.dev/nygenome/Lancet2/blob/db225350de4da2a125694e03e01dbc006a9865fc/website/static/img/chr14_72547800-72548098_c1_after_compression.png)

The fourth image below highlights all the assembly path flows taken through the graph.
![path_flow](https://github.dev/nygenome/Lancet2/blob/db225350de4da2a125694e03e01dbc006a9865fc/website/static/img/chr14_72547800-72548098_c1_path_flow.png)
