# Get the unbiased estimate of the correlation matrix for the
# population sample of the 21 Vulpes vulpes specimens

R.unb <- cor( vulpes21x29.matrix )

# Save the R.unb
devtools::use_data( R.unb )
