Model of   [Craig et al  *Bull Math Biol* **78** 2304-2357 (2016)](https://www.ncbi.nlm.nih.gov/pubmed/27324993) in the R package myelo accessed via shinyapp.io 

Move the *Move the Chemo Dose* slider and notice how the `renderPlot` expression is automatically re-evaluated when its dependant, `input$mult`, changes, causing a new plot of marrow cell state dynamocs to be rendered.
