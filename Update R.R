#This code automatically updates the R software installed in the computer
#It iinstalls the latest version of R and transfers the packages into the new version
if(!require(installr)) { 
install.packages("installr"); require(installr)} 
updateR()