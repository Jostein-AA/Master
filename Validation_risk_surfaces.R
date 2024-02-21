#Library needed to make GIFs
library(magick)

dir = "./Plots/scenario_1/"
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

image_write(my.animation, "test3.gif")
