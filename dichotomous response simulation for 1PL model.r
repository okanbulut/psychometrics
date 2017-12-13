########################################################
#IRT dichotomous response simulation                   #
#the function prob.irt takes the following argumensts: #
#abilities = vector of abilities.                      #
#b = vector of difficulties.                           #
########################################################

prob.irt = function(abilities,b,le) {

#makes abilities a column vector and the b parameters a row vector.
abilities = matrix(abilities,,1)
b = matrix(b,1,)
le = matrix(le,1,)

#obtains the number of examinees by looking at the length of the ability vector.
number.examinees = length(abilities)

#obtains the number of items by looking at the lenght of b vector.
number.items = length(b)

#calculates probabiities for the 1 PL.
      {
     #calculate the logit for the 1 PL: theta - b  (in matrix form)
    logit = matrix(abilities,number.examinees,number.items) -
            matrix(b,number.examinees,number.items,byrow = T)+ matrix(le,number.examinees,number.items,byrow = T)
    #convert logits to odds.
    odds = exp(logit)
    #convert odds to probabilities.
    prob1PL = odds/(1 + odds)
    #outputs matrix of probabilities for the 1 PL.
    return(prob1PL) }
      #close prob.irt function.
      }
#=====================================

##============================================================================
#STEP 4: CALCULATE THE PROBABILITY OF CORRECT RESPONSE.
#read population parameters from step 1 saved files.

population.difficulties = read.table(file="C:/R code/population_item_difficulty.dat",header=FALSE,sep="")
item.difficulties = population.difficulties[,2]

#apply testlet function and obtain the testlet design vector.
testlet1 = testlet(number.items=60, number.testlets=10, testlet.size=3, testlet.effect=0.5)

#create an item.testlet matrix.
item.testlet <-cbind(item.difficulties,testlet1)

#read abilities from population ability matrix.
population.abilities = read.table(file="C:/R code/population_abilities.dat",header=FALSE,sep="")
abilities = population.abilities[,2]


#======================================================================
### STEP 5 = OBTAIN 0/1 ITEM RESPONSES GIVEN PROBABILITIES.

iteration = 100        #define iteration.

for (k in 1:iteration)   {

            #uses the function prob.irt to obtain probabilities for the 1 PL. 20 different set of probabilities.
            probabilities = prob.irt(abilities,b = item.testlet[,1], le = item.testlet[,2])

            #generate a random uniform numbers between 0 and 1 for each response
            #to be used as a cut-off for the 0/1 response.
            #we have to obtain as many numbers as there are probabilities.
            numbers = runif(length(probabilities),min=0,max =1)

            #converts the numbers into a matrix with dimensions examinees X items.
            number.examinees = dim(probabilities)[1]      #row number.
            number.items = dim(probabilities)[2]          #column number.
            numbers = matrix(numbers,number.examinees,number.items)

           #compares numbers with probabilities.
           #if probability >= number then TRUE, if probability < number then FALSE.
           test = (probabilities >= numbers)

           #convert TRUE/FALSE logical values into 0/1 responses (multiplying by 1 does the trick).
           responses = test*1

           #save responses in a format that can be read by BILOG.
           #it is necessary to add a respondent identification number to the responses
           #for the analysis with BILOG.
           #the argument sp = "\t" saves as a tab delimited file.
           #examinee.id = 1:dim(responses)[1]

           #save simulated response data.
           filename = paste("testlet.condition2","iteration",k,".dat",sep="")
           write.table (responses,file = filename,sep = " ",row.names=FALSE,col.names=FALSE)

           #the end of the population parameter input.
          }

##############################################################################################
