
library(shinythemes);

shinyUI(
    fluidPage(theme = shinytheme("cosmo"),
              includeScript('www/tools.js'),
              ##css
              tags$head(tags$title("StressTest"),
                        tags$link(rel = "stylesheet",
                                  type = "text/css",
                                  href = "styles.css"),
                        tags$link(rel = "stylesheet",
                                  type = "text/css",
                                  href = "//fonts.googleapis.com/css?family=Oswald"),
                        tags$link(rel = "stylesheet",
                                  type = "text/css",
                                  href = "//fonts.googleapis.com/css?family=Lora"),
                        tags$script(src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML",
                                    type="text/javascript")
                        ),

              ##title box
              withTags({
                  div(class="cheader",
                      "Clinical Trial Design Robustness Stress Test",
                      tags$button(
                               id = 'close',
                               type = "button",
                               class = "btn action-button",
                               onclick = "setTimeout(function(){window.close();},500);",
                               "Exit",
                               style="float: right;
                                      background-image: url(texturebg2.jpg);"
                           ))
              }),

              ##main page
              uiOutput("mainpage"),

              ##foot
              withTags({
                  div(class="cfooter", "OLSSOL Project")})
              )
)
