import random, math, time

#Sigmoid Function
#Default: calculates an output between 0 and 1 from any net input
#Derivative: calculates the gradient at a given output value
def sigmoid(x, deriv=False):
    if deriv:
        return x * (1 - x)
    return 1 / (1 + pow(math.e, -x))

#Mean Squared Error function
#Calculates the MSE between the target and actual outputs
def squaredError(target, output):
        return 0.5 * math.pow(target-output, 2)

#Neuron Class
class Neuron:
    #Creates a weight for each neuron in the previous layer plus a bias weight
    def __init__(self, n_weights):
        self.weights = [random.uniform(-1,1) for i in range(n_weights + 1)]
        self.weight_gradients = [0.0 for i in range(n_weights + 1)]
    #Calculates the net input to the neuron from the inputs and weights
    #Then it calculates the neuron's output using the sigmoid function
    def computeOutput(self, inputs):
        self.net_input = self.weights[-1] + sum([inputs[i]*self.weights[i] for i in range(len(inputs))])
        self.output = sigmoid(self.net_input)
        return self.output

#Layer Class
class Layer:
    #Creates the layers of neurons with the correct number of weights
    def __init__(self, n_neurons, weights_per_neuron):
        self.n_neurons = n_neurons
        self.weights_per_neuron = weights_per_neuron
        self.neurons = [Neuron(weights_per_neuron) for i in range(n_neurons)]

#Network Class
class Network:
    #Creates a network from the structure array (array includes input neurons)
    def __init__(self, structure):
        self.structure = structure
        self.layers = [Layer(structure[i], structure[i-1]) for i in range(1,len(structure))]
    #computes the outputs of each neuron
    def feedforward(self, inputs):
        self.layer_outputs = [inputs]
        for l in self.layers:
            self.layer_outputs.append([n.computeOutput(self.layer_outputs[-1]) for n in l.neurons])
        return self.layer_outputs[-1]
    def computeWeightChanges(self, targets):
        for l in range(len(self.layers)):
            l = -(l+1)
            for n in range(len(self.layers[l].neurons)):
                if l == -1:
                    self.layers[l].neurons[n].node_delta = -(targets[n] - self.layer_outputs[l][n]) * sigmoid(self.layer_outputs[l][n], deriv=True)
                else:
                    self.layers[l].neurons[n].node_delta = sum([(self.layers[l+1].neurons[j].node_delta * self.layers[l+1].neurons[j].weights[n]) for j in range(len(self.layers[l+1].neurons))]) * sigmoid(self.layer_outputs[l][n], deriv=True)
                for w in range(len(self.layers[l].neurons[n].weights)-1):
                    self.layers[l].neurons[n].weight_gradients[w] = self.layers[l].neurons[n].node_delta * self.layer_outputs[l-1][w]
                self.layers[l].neurons[n].weight_gradients[-1] = self.layers[l].neurons[n].node_delta * 1.0
    def changeWeights(self, learning_rate):
        for l in range(len(self.layers)):
            for n in range(len(self.layers[l].neurons)):                
                for w in range(len(self.layers[l].neurons[n].weights)):
                    self.layers[l].neurons[n].weights[w] -= learning_rate * self.layers[l].neurons[n].weight_gradients[w]
    def train(self, input_dataset, target_dataset, learning_rate, cycles=1, display_interval=1.0):
        time_start = time.time()
        time_then = time_start
        for cycle in range(cycles):
            for case in range(len(input_dataset)):
                #TIMING CODE BEGINNING
                # note that the timing code takes up ~0.0002s per cycle for a [784,30,10] network
                time_now = time.time()
                if time_now - time_then > display_interval:
                    time_then = time_now
                    print("Cycle "+str(cycle + 1)+"/"+str(cycles)+"    Training case "+str(case + 1)+"/"+str(len(input_dataset)))
                #TIMING CODE END
                self.feedforward(input_dataset[case])
                self.computeWeightChanges(target_dataset[case])
                self.changeWeights(learning_rate)
        time_taken = time.time() - time_start
        print("Training",str(cycles),"cycles of",str(len(input_dataset)),"cases complete in",str(round(time_taken, 5))+"s.")
        return time_taken
                    
                
        

#setup and operations:
random.seed(1)
net = Network([2,2,2,2])
print([[[w for w in n.weights] for n in l.neurons] for l in net.layers])

print()
net.feedforward([0.1,0.9])
print(net.layer_outputs)

print()
net.computeWeightChanges([0.05,0.10])
net.changeWeights(0.5)

net.train([[0.1,0.9],[0.9,0.1]], [[0.05,0.10],[0.10,0.05]], 0.5, cycles=50000)

#a function for easily showing the value of a certain property for each neuron in the network
def show(prop):
    global net
    for l in net.layers:
        for n in l.neurons:
            print(eval("n."+prop))
