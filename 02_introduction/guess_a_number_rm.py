# This is a guess the number game.
import random

guessesTaken = 0

print('Welcome to my Game')
print('I am thinking of a number between 0 and 100!')

number = random.randint(0, 100)


while guessesTaken < 4:
    numTries = 3 - guessesTaken
    numTriesStr = str(numTries)
    print('you have '+ numTriesStr +' tries')
    print('Take a guess.')
    guess = input()
    guess = int(guess)

    guessesTaken = guessesTaken + 1
    if guess < number:
        print('Your guess is too low.')

    if guess > number:
        print('Your guess is too high.')

    if guess == number:
        break

if guess == number:
        guessesTaken = str(guessesTaken)
        print('Good job! You guessed my number in ' + guessesTaken + ' guesses!')

if guess != number:
        number = str(number)
        print('Nope. The number I was thinking of was ' + number)
