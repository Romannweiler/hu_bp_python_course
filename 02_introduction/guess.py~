from random import randint
from sys import exit

if __name__ == "__main__":
	retry=True
	while retry:
		z=randint(0,100)
		print z
		#userinput
		guesscounter=0
		input2=-2
		while input2!=z:
			input2=raw_input('Guess my number!')
			#verwertbarer input?
			neueingeben=True
			while neueingeben:
				try:
					input2=int(input2)
					neueingeben=False
				except:
					input2=raw_input('Your input was not useable. Reinput it please. This will not count towards your total count')
			guesscounter+=1
			if input2>z:
				print 'My number is lesser than yours.'
			elif input2<z:
				print 'My number is greater than yours.'
			elif input2==z:
				print 'You guessed my number!'
				print 'You needed '+str(guesscounter)+' guess(es).'
				break
		#wanna try again?
		input1='0'
		while ( not isinstance(input1 ,int)) and (input1 !=0 or input1!=1):
			input1=raw_input('wanna play again? 1=yes,0=no')
			try:
				input1=int(input1)
			except:
				input1='0'
		if input1==False:
			retry = False
	exit()
