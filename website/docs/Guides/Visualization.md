# Visualization

By adding the ```--graphs-dir``` tag and providing a directory where you want the graphs to be placed, Lancet will construct deBruijn graphs for each window inspected and place them into the directory as .dot files. NOTE: whichever directory is given as the graphs-dir will be cleared so be mindful of what directory you provide.

```bash
./lancet2 pipeline -t tumor.bam -n normal.bam -r ref.fasta --graphs-dir ./graph_dir
```
The above command will export the DeBruijn graph after every stage of the assembly (low covergae removal, tips removal, compression) to the following set of files:

1. chr:start-end_cX_before_pruning.dot
2. chr:start-end_cX_after_pruning.dot
3. chr:start-end_cX_path_flow.dot

Where X is the number of the correspending connected component (in most cases only one). 
These files can be rendered using the utilities available in the [Graphviz](http://www.graphviz.org/) visualization software package. Specifically we recommend using the **sfdp** utlity which draws undirected graphs using the ``spring'' model and it uses a multi-scale approach to produce layouts of large graphs in a reasonably short time.

```
sfdp -Tpdf example_file.dot -O
```

The above command will create a example_file.dot.pdf file that shows the graph. For large graphs, Adobe Acrobat Reader may have troubles rendering the graph in which case we recommend opening the PDF file using the "Preview" image viewer software available in MacOS.

# PROVIDE EXAMPLE IMAGES
What is a good window to see similar features as was in the first lancet pictures
