# Java skeleton for duck hunt dd2380

## Compile
```
javac *.java
```
## Run
The agent can be run in two different modes:
2. Server - act as the judge by sending predefined observations one at a time
  and asking the client to respond 
3. Client - get observations from standard input and output actions to
  standard output (this is the default mode)

The server and client can be run in separate terminals and communicate
through pipes. Create the pipes first. (we recommend cygwin for windows users)
mkfifo player2server server2player

Terminal 1:
```
java Main verbose server < player2server > server2player
```

Terminal 2:
```
java Main verbose > player2server < server2player
```

Or you may run both instances in the same terminal
```
java Main server < player2server | java Main verbose > player2server
```

You can test a different environment like this
```
java Main server load ParadiseEmissions.in < player2server | java Main verbose > player2server
```

## How it works:
### Guessing:
* first round is random
* other rounds calculate the model that gives the max likelihood of the unknown bird, with 2 constraints:
  * logProb for that species must be over a threshold (ex [-200 -1300 -1400], choose 0)
  * logProb must be significantly higher than the others (ex [-200 -205 -1400], do not choose 0 because 1 is too near)
* if the bird is unknown or has a logProb of being a black stork that's higher than the black stork threshold add that bird to a list of risky shots
* for every unknown bird consider guessing randomly with an unseen species, but no more than 2 random guesses per unseen species (ex [-1 2 2 -1 -1 -1 -1 3 4] with unseen species 0 and 5 -> [0 2 2 5 0 5 -1 3 4])
* last round no random guessing
### Shooting
* use the species model to identify a bird
* use the species model to predict next move (skip calulating for black storks)
* choose the most probable move among all the most probable moves of the birds
* shoot if prob is greater than a threshold
### Reval
* for every revealed bird add its sequence to the list of sequences of its species
* retrain the model of every species using all the sequences
### Other things
* smothing in the alpha pass to avoid dividing by zero
* 3 states are enough because not all species exhibit all 5 behaviors
* 100 iterations max for Baum Welch are enough
* the models of the birds created during the rounds are useless, it's better to use the model of the species to shoot
