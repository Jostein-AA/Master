#Library needed to make GIFs
library(magick)

make_GIF <- function(dir, gif_name){
  files <- list.files(path = dir, pattern = "*.png")
  for(i in 1:length(files)){
    if(paste(toString(i), ".png", sep = "") %in% files){
      if(i == 1){
        img = c(image_read(paste(dir, 
                                 paste(toString(i), ".png", sep = ""), 
                                 sep = "")))
      } else{
        img = c(img, image_read(paste(dir, 
                                      paste(toString(i), ".png", sep = ""), 
                                      sep = "")))
      }
    } else{
      print(paste("No such file:", paste(toString(i), ".png", sep = "")))
    }
  }
  
  image_append(image_scale(img, "x200"))
  
  my.animation <- image_animate(image_scale(img, 
                                            "400x400"),
                                fps = 1,
                                dispose = "previous")
  print("Writing")
  image_write(my.animation, gif_name)
  
}

dir = "./Plots/scenario_1/"
gif_name = "sc1.gif"

make_GIF(dir, gif_name)

dir = "./Plots/scenario_3/"
gif_name = "sc3.gif"

make_GIF(dir, gif_name)

dir = "./Plots/scenario_5/"
gif_name = "sc5.gif"

make_GIF(dir, gif_name)

dir = "./Plots/scenario_7/"
gif_name = "sc7.gif"

make_GIF(dir, gif_name)

dir = "./Plots/scenario_9/"
gif_name = "sc9.gif"

make_GIF(dir, gif_name)

dir = "./Plots/scenario_11/"
gif_name = "sc11.gif"

make_GIF(dir, gif_name)



