# factoring
 A suite of factoring algorithms I decided to make, for Terra's week 4 "Splitting" theme.

For the theme of splitting, because my special interest was "factoring large integers" among other things, I decided to actually get to implementing a factoring suite in python (yes yes I know, python is slow, yada yada yada, but JS has no easy arbitrary precision, and I don't know how to write C-family.) I originally wasn't going to include so many algorithms, but changed my mind and added basically every one of them below quadratic sieve.

The technical difficulty was immense in this project. Much of the project above implementational details were heavily number theory based and thus a lot of unlogged time was spent on researching, implementing, cross-checking against the papers, etc. Unfortunately that's not really loggable time, so whatever. 

This project was somewhat rushed due to exam week being on the week the project was made during. Hackclub Haven also fell during the weekend in which I was an org, which made me effectively unable to code for the majority of daytime during those two days. I wanted to include SO MANY other algorithms (like Dixon's, William's P+1, Hart's one-line, older and curiosity algorithms like CFRAC and SQUFOF, etc...) but it just wasn't really enough time...

The program will output `C!<number>` if the last cofactor is known composite!

# Running

TO RUN:

- Download the file `factor-linux`
- Run `./factor-linux` (`chmod 777 ./factor-linux` if any permission errors occur)

# AI use

AI was used to explain pollard's P-1, ECM, and MPQS algorithms to me. 
AI was used to DEBUG on a number-theory level after I failed to find the issue after bashing my head against the code for 10-15 minutes.
AI was NOT used to write any direct code.

# Reviewers

For reviewers: here are some sample numbers to feed:
- 19999999999999999999999999999999999111111111
- 117634454544619836473649514876284198472975319873
- 1384756789238574696814975978659821658727059786598476283676729176458197681
- 3250909524224980026658272626172851784216727881201572167748724186378794282641628002574416019948169976992893495329618433676502898369445077616884991205431