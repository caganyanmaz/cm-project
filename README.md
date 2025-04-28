## The Continuous Maths Vacation Project 

I used Newton's Method with backtracking (applying Armijo Rule) to solve the given optimization problem.

# How to run the program 

Run make.

If you don't have make installed, you can run the following code:

```
    gcc main.c matrix.c -lm -o main
```

This will create a executable called main in the file directory. 

In order to use the program, you need to include the parameters of the images in a file named "cm-sheet5-X.txt" and the labels of the images in a file named "cm-sheet5-y.txt" in the project directory. If you want to change these names, you need to modify it in the source code.

After running the program, it'll calculate weights and stop when it reaches one of the end conditions (good enough estimate, Armijo Rule fails, too many steps), it'll print the accuracy and save the weights in a file called "cm-sheet5-w.txt" in project directory. If you call the program again, it'll load weights of this file for further optimization. If you don't want that to happen, you need to rename the file.

There's also a bug where when I first run it, it ends the optimization process early and get a worse accuracy. But if you run the program a multiple times it gets better, although a terminal condition was reached in the first run. I don't know why that happens. If you run it 3-4 times, you should get ~95% accuracy

The weights in my machine have accuracy 94.8%, but it should vary as initially the program assigns weights randomly using a Normal Distribution estimate.


