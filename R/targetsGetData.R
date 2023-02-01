tar_option_set(packages = c("tidyverse", "lubridate"))

target_getData <- list(
  tar_target(d0_file, "./dataIn/JFIP_ElectrofishingData_Converted.csv", format = "file"),
  tar_target(d0, {
    read_csv(d0_file)
  }),
  
  tar_target(d0_2021_file, "./dataIn/JFIP_ElectrofishingData_2021_ConvertedCoordinates.csv", format = "file"),
  tar_target(d0_2021, {
    read_csv(d0_2021_file)
  })
  
  # tar_target(sampleDateMedian_file, "./dataOut/sampleDateMedian.csv", format = "file"),
  # tar_target(target_sampleDateMedian, {
  #   read_csv(sampleDateMedian_file)
  # })
)
  # create_output <- function(file) {
#   data <- read.csv(file)
#   output <- head(data)
#   write.csv(output, paste0("./dataForModels/", fileName, ".csv"))
#   "output.csv"
# }

# list(
#   tar_target(name = input_d0, command = "./dataIn/JFIP_ElectrofishingData_Converted.csv", format = "file"),
#   tar_target(name = output, command = create_output(input_d0, "d0In"), format = "file"),
#   
#   tar_target(name = input_d02021, command = "./dataIn/JFIP_ElectrofishingData_2021_ConvertedCoordinates.csv", format = "file"),
#   tar_target(name = output, command = create_output(input_d02021, "d0In2021"), format = "file")
# )
