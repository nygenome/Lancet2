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

Below is an example of what the generated graphs may look like. The first image is the raw graph before removing the low coverage nodes, the second image is before graph compression, the third image is after graph compression, and the fourth image highlights all the path flows taken through the graph during assembly. The blue nodes are k-mers shared by both tumor and normal; the white nodes are k-mer with low support (likely sequencing errors); the green nodes are k-mers only present in the normal; the red nodes are k-mers only present in the tumor.:

![raw_graph](https://github.com/nygenome/Lancet2/tree/main/website/static/img/chr14_72547800-72548098_c0_raw_graph.png)
![before_compression](https://github.com/nygenome/Lancet2/tree/main/website/static/img/chr14_72547800-72548098_c1_before_compression.png)
![after_compression](https://github.com/nygenome/Lancet2/tree/main/website/static/img/chr14_72547800-72548098_c1_after_compression.png)
![path_flow](https://github.com/nygenome/Lancet2/tree/main/website/static/img/chr14_72547800-72548098_c1_path_flow.png)
