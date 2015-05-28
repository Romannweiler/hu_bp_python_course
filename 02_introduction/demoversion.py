from random import randint
from sys import exit


def solstrat(z,a,b,count):

	ledge=a
	redge=b
	randguess=randint(a,b)
	print 'The player guesses '+str(randguess)+'.'	
	counter=count+1
	if randguess<z:
		print 'The guessed number is lesser than the solution.'
		ledge=randguess
		solstrat(z,ledge,redge,counter)
	elif randguess>z:
		print 'The guessed number is greater than the solution.'
		redge=randguess
		solstrat(z,ledge,redge,counter)
	elif randguess==z:
		print 'The player guessed the number!'
		print 'That number was '+str(randguess)
		print 'The Player needed '+str(counter)+' guess(es).'
	return None



if __name__ == "__main__":
	print 'Welcome to the demoversion of "Guess a number!".\n For demonstration purposes the computer will play against itself.'
	z=randint(0,100)
	print 'First of all, the programm guesses a number, which will be the solution to this game.\n The number is: '+str(z)
	print 'Now the player will start:'
	guesscounter=0
	ledge=0
	redge=100
	solstrat(z,ledge,redge,guesscounter)
	
	print'The game was finished successfully. For the full version, see guess.py. Goodbye.'

	exit()
