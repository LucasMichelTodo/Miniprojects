library(tsne)
library(plotly)
Sys.setenv("plotly_username"="LucasMichelTodo")
Sys.setenv("plotly_api_key"="P7yXFaubnGsAb8WwemL7")


trans <- read.csv2("/media/lucas/Disc4T/Projects/PhD_Project/External_Data/3D7variantome_10g12b3d7b_trans.csv", sep = ",")

trans <- trans[,c(5,2,3,4)]

## 3D Plot
p <- plot_ly(trans, x = ~X1.2B, y = ~X10G, z = ~X3D7.B) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = '1.2B'),
                      yaxis = list(title = '10G'),
                      zaxis = list(title = '3D7.B')))

# Create a shareable link to your chart
# Set up API credentials: https://plot.ly/r/getting-started
chart_link = api_create(p, filename="scatter3d-basic")
chart_link

