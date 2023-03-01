if (!require(styler)) {
  install.packages("styler")
}

print("Tidying code...")
style_dir(exclude_files = "sandbox.r")
