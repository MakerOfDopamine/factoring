# factoring
 A suite of factoring algorithms I decided to make, for Terra's week 4 "Splitting" theme.

For the theme of splitting, because my special interest was "factoring large integers" among other things, I decided to actually get to implementing a factoring suite in python (yes yes I know, python is slow, yada yada yada, but JS has no easy arbitrary precision, and I don't know how to write C-family.) I originally wasn't going to include so many algorithms, but changed my mind and added basically every one of them below quadratic sieve.

The technical difficulty was immense in this project. Much of the project above implementational details were heavily number theory based and thus a lot of unlogged time was spent on researching, implementing, cross-checking against the papers, etc. Unfortunately that's not really loggable time, so whatever. 

This project was somewhat rushed due to exam week being on the week the project was made during. Hackclub Haven also fell during the weekend in which I was an org, which made me effectively unable to code for the majority of daytime during those two days. I wanted to include SO MANY other algorithms (like Dixon's, William's P+1, Hart's one-line, older and curiosity algorithms like CFRAC and SQUFOF, etc...) but it just wasn't really enough time...

# Running

TO RUN:

- Download the source code
- Unzip the source code
- and run "python3 factor.py" from within the directory.
- Numbers to try are in the bottom of the file.

# AI use

AI was used to explain pollard's P-1, ECM, and MPQS algorithms to me. 
AI was used to DEBUG on a number-theory level after I failed to find the issue after bashing my head against the code for 10-15 minutes.
AI was NOT used to write any direct code.