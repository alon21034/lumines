/* (function() { // use closure technique to hide our functions */ 
/*******************************
 * We are using NEAT algorithm.
 * This script is based on http://pastebin.com/ZZmSNaHX, which is writen in lua.
 *******************************/

  function sigmoid(x) {
    return 2 / (1 + Math.exp(-4.9*x)) - 1;
  }

  var Population = 300;
  var DeltaDisjoint = 2.0;
  var DeltaWeights = 0.4;
  var DeltaThreshold = 1.0;

  var StaleSpecies = 15;

  var MutateConnectionsChance = 0.25;
  var PerturbChance = 0.90;

  var CrossoverChance = 0.75;
  var LinkMutationChance = 2.0;
  var NodeMutationChance = 0.50;
  var BiasMutationChance = 0.40;
  var StepSize = 0.1;
  var DisableMutationChance = 0.4;
  var EnableMutationChance = 0.2;

  var TimeoutConstant = 50;

  var MaxNodes = 10000000;

  var inputKeys = ['left', 'up', 'right', 'down', 'space'];

  var gridSize = 16 * 10;
  var bricks = 4 * (5 + 1); // one current brick, 5 waiting list

  var nOutput = inputKeys.length;
  // for each grid, save color and darkness,
  // 2 is for x, y position of current brick,
  // 1 + 1 is for bar.x and constant bias
  var nInput = gridSize * 2 + bricks + 2 + 1 + 1;

  var NUM_COL = 16;
  var NUM_ROW = 10;

  var readInputs = function(game, map, brick, bar) {
    var inputs = [];

    // 1. read value of grids
    for (var col = 0; col < NUM_COL; ++ col) {
      for (var row = 0; row < NUM_ROW; ++ row) {
        switch (map.grid[col][row]) {
          case 0:
            inputs.push(0);
            inputs.push(0);
            break;
          case 1:
            inputs.push(-1);
            inputs.push(0);
            break;
          case 2:
            inputs.push(1);
            inputs.push(0);
            break;
          case 3:
            inputs.push(-1);
            inputs.push(1);
            break;
          case 4:
            inputs.push(1);
            inputs.push(1);
            break;
        }
      }
    }

    // 2. read bricks
    for (var i = 0; i < brick.color.length; ++ i) {
      inputs.push(brick.color[i] * 2 - 1);
    }
    for (var i = 0; i < brick.nextFive.length; ++ i) {
      for (var j = 0; j < brick.nextFive[i].length; ++ j) {
        inputs.push(brick.nextFive[i][j] * 2 - 1);
      }
    }
    /*
    inputs = inputs.concat(brick.color);
    for (var i = 0; i < 5; ++ i) {
      inputs = inputs.concat(brick.nextFive[i]);
    }
    */

    inputs.push(brick.x);
    inputs.push(brick.y);

    inputs.push(bar.x);

    if (inputs.length != nInput - 1) {
      console.log("WARNING! inputs.length does not math (in readInputs)");
    }

    return inputs;
  }

  /**
   * define NEAT functions
   */

  /**
   * a global counter, 
   *  innovation.next() increases the counter by 1
   *  innovation.value() returns current value
   */
  var Innovation = function(v) {
    var value = v;
    this.next = function() {
      value = value + 1;
      return value;
    }
    this.value = function() {
      return value;
    }
  }

  var innovation = new Innovation(nOutput);

  var Pool = function() {
    this.species = [];
    this.generation = 0;
    innovation = new Innovation(nOutput);
    this.currentSpecies = 0;
    this.currentGenome = 0;
    this.currentFrame = 0;
    this.maxFitness = 0;
  };

  Pool.prototype.initialize = function() {
    this.species = [];
    this.generation = 0;
    innovation = new Innovation(nOutput);
    this.currentSpecies = 0;
    this.currentGenome = 0;
    this.currentFrame = 0;
    this.maxFitness = 0;
    for (var i = 0; i < Population; ++ i) {
      var basic = Genome.getBasicGenome();
      this.addToSpecies(basic);
    }
  }

  Pool.prototype.rankGlobally = function() {
    var global = [];
    for (var i = 0; i < this.species.length; ++ i) {
      var sp = this.species[i];
      for (var j = 0; j < sp.genomes.length; ++ j) {
        global.push(sp.genomes[j]);
      }
    }
    global.sort(function(a, b) { return a.fitness - b.fitness; });
    for (var i = 0; i < global.length; ++ i) {
      global[i].globalRank = i;
    }
  }

  Pool.prototype.totalAverageFitness = function() {
    var total = 0;
    for (var i = 0; i < this.species.length; ++ i) {
      var sp = this.species[i];
      total += sp.averageFitness;
    }
    return total;
  }

  Pool.prototype.cullSpecies = function(cutToOne) {
    for (var i = 0; i < this.species.length; ++ i) {
      var sp = this.species[i];
      sp.genomes.sort(function(a, b) { return b.fitness - a.fitness; });

      var remaining = Math.ceil(sp.genomes.length / 2);
      if (cutToOne) remaining = 1;
      while (sp.genomes.length > remaining) {
        sp.genomes.pop();
      }
    }
  }

  Pool.prototype.removeStaleSpecies = function() {
    var survived = [];
    for (var i = 0; i < this.species.length; ++ i) {
      var sp = this.species[i];
      sp.genomes.sort(function (a, b) {return b.fitness - a.fitness;});

      if (sp.genomes[0].fitness > sp.topFitness) {
        sp.topFitness = sp.genomes[0].fitness;
        sp.staleness = 0;
      } else {
        sp.staleness ++;
      }

      if (sp.staleness < StaleSpecies || sp.topFitness > this.maxFitness) {
        survived.push(sp);
      }
    }
    this.species = survived;
  }

  Pool.prototype.removeWeakSpecies = function() {
    var survived = [];
    var sum = this.totalAverageFitness();
    for (var i = 0; i < this.species.length; ++ i) {
      var sp = this.species[i];
      var breed = Math.floor(sp.averageFitness / sum * Population);
      if (breed >= 1) {
        survived.push(sp);
      }
    }
    this.species = survived;
  }

  Pool.prototype.addToSpecies = function(child) {
    var foundSpecies = false;
    for (var i = 0; i < this.species.length; ++ i) {
      var sp = this.species[i];
      if (!foundSpecies && Genome.sameSpecies(child, sp.genomes[0])) {
        sp.genomes.push(child);
        foundSpecies = true;
        break;
      }
    }
    if (!foundSpecies) {
      var sp = new Species();
      sp.genomes.push(child);
      this.species.push(sp);
    }
  }

  Pool.prototype.newGeneration = function() {
    this.cullSpecies(false); // cull the buttom half of each species
    this.rankGlobally();
    this.removeStaleSpecies();
    this.rankGlobally();
    for (var i = 0; i < this.species.length; ++ i) {
      var sp = this.species[i];
      sp.calculateAverageFitness();
    }
    this.removeWeakSpecies();
    var sum = this.totalAverageFitness();
    var children = [];
    for (var i = 0; i < this.species.length; ++ i) {
      var sp = this.species[i];
      var breed = Math.floor(sp.averageFitness / sum * Population) - 1;
      for (var j = 0; j < breed; ++ j) {
        children.push(sp.breedChild());
      }
    }
    this.cullSpecies(true);
    while (children.length + this.species.length < Population) {
      var i = Math.floor(Math.random() * this.species.length);
      var sp = this.species[i];
      children.push(sp.breedChild());
    }
    for (var i = 0; i < children.length; ++ i) {
      var child = children[i];
      this.addToSpecies(child);
    }
    this.generation ++;
    // TODO save this generation

    this.save("backup");
  }

  Pool.prototype.save = function(filename) {
    var storage = {};

    storage.generation = this.generation;
    storage.maxFitness = this.maxFitness;
    storage.species = [];
    for (var i = 0; i < this.species.length; ++ i) {
      storage.species.push(JSON.parse(this.species[i].pack()));
    }

    // localStorage.setItem(filename, JSON.stringify(storage));
    neatDataBase.saveGeneration(this);
  }

  Pool.prototype.load = function(id) {
    var caller = this;
    neatDataBase.getGeneration(id, function(storage) {
      caller.generation = storage.generation;
      caller.maxFitness = storage.maxFitness;
      caller.species = [];
      for (var i = 0; i < storage.species.length; ++ i) {
        caller.species.push(Species.unpack(JSON.stringify(storage.species[i])));
      }
      console.log("got it!");
    });
    /*
    var storage = JSON.parse(localStorage.getItem(filename));

    this.generation = storage.generation;
    this.maxFitness = storage.maxFitness;
    this.species = [];
    for (var i = 0; i < storage.species.length; ++ i) {
      this.species.push(Species.unpack(JSON.stringify(storage.species[i])));
    }
    */
  }

  var Species = function() {
    this.topFitness = 0;
    this.staleness = 0;
    this.genomes = [];
    this.averageFitness = 0;
  }

  Species.prototype.pack = function() {
    var storage = {};
    storage.topFitness = this.topFitness;
    storage.staleness = this.staleness;
    storage.genomes = [];
    for (var i = 0; i < this.genomes.length; ++ i) {
      storage.genomes.push(JSON.parse(this.genomes[i].pack()));
    }
    return JSON.stringify(storage);
  }

  Species.unpack = function(packed) {
    var storage = JSON.parse(packed);
    var sp = new Species();
    sp.topFitness = storage.topFitness;
    sp.staleness = storage.staleness;
    sp.genomes = [];
    for (var i = 0; i < storage.genomes.length; ++ i) {
      sp.genomes.push(Genome.unpack(JSON.stringify(storage.genomes[i])));
    }
    return sp;
  }

  Species.prototype.calculateAverageFitness = function() {
    var total = 0;
    for (var i = 0; i < this.genomes.length; ++ i) {
      var genome = this.genomes[i];
      total += genome.globalRank;
    }
    this.averageFitness = total / (this.genomes.length);
  }

  Species.prototype.breedChild = function() {
    var child = [];
    var nGenomes = this.genomes.length;
    if (Math.random() < CrossoverChance) {
      var g1 = this.genomes[Math.floor(Math.random() * nGenomes)];
      var g2 = this.genomes[Math.floor(Math.random() * nGenomes)];
      child = Genome.crossover(g1, g2);
    } else {
      var g = this.genomes[Math.floor(Math.random() * nGenomes)];
      child = g.clone();
    }
    child.mutate();
    return child;
  }

  var Genome = function() {
    this.genes = [];
    this.fitness = 0;
    this.adjustedFitness = 0;
    this.network = {};
    this.maxNeuron = 0;
    this.globalRank = 0;
    this.mutationRates = {};

    this.mutationRates["connections"] = MutateConnectionsChance;
    this.mutationRates["link"] = LinkMutationChance;
    this.mutationRates["bias"] = BiasMutationChance;
    this.mutationRates["node"] = NodeMutationChance;
    this.mutationRates["enable"] = EnableMutationChance;
    this.mutationRates["disable"] = DisableMutationChance;
    this.mutationRates["step"] = StepSize;
  }

  Genome.prototype.display = function() {
    var network = this.network;
  }

  Genome.prototype.pack = function() {
    var storage = {};
    storage.fitness = this.fitness;
    storage.maxNeuron = this.maxNeuron;
    storage.mutationRates = this.mutationRates;
    storage.genes = [];
    for (var i = 0; i < this.genes.length; ++ i) {
      storage.genes.push(JSON.parse(this.genes[i].pack()));
    }
    return JSON.stringify(storage);
  }

  Genome.unpack = function(packed) {
    var g = new Genome();
    var storage = JSON.parse(packed);
    g.fitness = storage.fitness;
    g.maxNeuron = storage.maxNeuron;
    g.genes = [];
    for (var i = 0; i < storage.genes.length; ++ i) {
      g.genes.push(Gene.unpack(JSON.stringify(storage.genes[i])));
    }
    return g;
  }

  Genome.prototype.generateNetwork = function() {
    this.network = new NeuronNetwork(this);
  }

  Genome.crossover = function(g1, g2) {
    if (g2.fitness > g1.fitness) {
      return Genome.crossover(g2, g1);
    }

    var child = new Genome();
    var innovations2 = [];
    for (var i = 0; i < g2.genes.length; ++ i) {
      var gene = g2.genes[i];
      innovations2[gene.innovation] = gene;
    }
    for (var i = 0; i < g1.genes.length; ++ i) {
      var gene1 = g1.genes[i];
      var gene2 = innovations2[gene1.innovation];
      if (gene2 !== undefined && gene2.enabled && Math.random() > 0.5) {
        child.genes.push(gene2.clone());
      } else {
        child.genes.push(gene1.clone());
      }
    }
    child.maxNeuron = Math.max(g1.maxNeuron, g2.maxNeuron);
    for (var mutation in g1.mutationRates) {
      child.mutationRates[mutation] = g1.mutationRates[mutation];
    }
    return child;
  }

  /**
   * randomly choose on neuron.
   * if nonInput == true, don't choose from input neurons.
   */
  Genome.prototype.randomNeuron = function(nonInput) {
    var candidates = [];
    if (!nonInput) {
      for (var i = 0; i < nInput; ++ i) {
        candidates[i] = true;
      }
    }
    for (var i = 0; i < nOutput; ++ i) {
      candidates[MaxNodes + i] = true;
    }
    for (var i = 0; i < this.genes.length; ++ i) {
      var gene = this.genes[i];
      if (!nonInput || gene.into >= nInput) {
        candidates[gene.into] = true;
      }
      if (!nonInput || gene.out >= nInput) {
        candidates[gene.out] = true;
      }
    }
    var count = 0;
    for (var idx in candidates) count ++;
    var n = Math.floor(Math.random() * count);
    for (var idx in candidates) {
      if (n == 0) {
        return idx;
      }
      n --;
    }
    return 0;
  }

  Genome.prototype.containsLink = function(link) {
    for (var i = 0; i < this.genes.length; ++ i) {
      var gene = this.genes[i];
      if (gene.into == link.into && gene.out == link.out) {
        return true;
      }
    }
    return false;
  }

  Genome.prototype.pointMutate = function () {
    var step = this.mutationRates["step"];
    for (var i = 0; i < this.genes.length; ++ i) {
      var gene = this.genes[i];
      if (Math.random() < PerturbChance) {
        gene.weight = gene.weight + Math.random() * step * 2 - step;
      } else {
        gene.weight = Math.random() * 4 - 2;
      }
    }
  }

  Genome.prototype.linkMutate = function(forceBias) {
    var neuron1 = this.randomNeuron(false);
    var neuron2 = this.randomNeuron(true);
    var link = new Gene();
    if (neuron1 < nInput && neuron2 < nInput) {
      // Both input nodes
      return;
    }
    if (neuron2 < nInput) {
      var t = neuron1;
      neuron1 = neuron2;
      neuron2 = t;
    }
    link.into = neuron1;
    link.out = neuron2;
    if (forceBias) {
      link.into = nInput - 1;
    }

    if (this.containsLink(link)) {
      return;
    }
    link.innovation = innovation.next();
    link.weight = Math.random() * 4 - 2;
    // console.log("new link: " + link.into + " " + link.out + " " + link.weight);
    this.genes.push(link);
  }

  /**
   * break a link into two links by adding an additional node
   * A --------> B becomes A --> C --> B
   */
  Genome.prototype.nodeMutate = function() {
    if (this.genes.length == 0) return; // ???

    this.maxNeuron ++;

    var i = Math.floor(Math.random() * this.genes.length);
    var gene = this.genes[i];
    if (!gene.enabled) return;

    var gene1 = gene.clone();
    gene1.out = this.maxNeuron - 1;
    gene1.weight = 1;
    gene1.innovation = innovation.next();
    gene1.enabled = true;
    this.genes.push(gene1);

    var gene2 = gene.clone();
    gene2.into = this.maxNeuron - 1;
    gene2.innovation = innovation.next();
    gene2.enabled = true;
    this.genes.push(gene2);
  }

  Genome.prototype.enableDisableMutation = function(enable) {
    var candidates = [];
    for (var idx in this.genes) {
      if (this.genes[idx].enabled === !enable) {
        candidates.push(this.genes[idx]);
      }
    }
    if (candidates.length == 0) return;
    var which = Math.floor(Math.random() * candidates.length);
    candidates[which].enabled = !candidates[which].enabled;
  }

  Genome.prototype.mutate = function() {
    for (var mutation in this.mutationRates) {
      var rate = this.mutationRates[mutation];
      if (Math.random() > 0.5) {
        this.mutationRates[mutation] = 0.95 * rate;
      } else {
        this.mutationRates[mutation] = 1.05263 * rate;
      }
    }

    if (Math.random() < this.mutationRates["connections"]) {
      this.pointMutate();
    }

    var p;
    p = this.mutationRates["link"];
    while (p > 0) {
      if (Math.random() < p) {
        this.linkMutate(false);
      }
      p -= 1;
    }

    p = this.mutationRates["bias"];
    while (p > 0) {
      if (Math.random() < p) {
        this.linkMutate(true);
      }
      p -= 1;
    }
    p = this.mutationRates["node"];
    while (p > 0) {
      if (Math.random() < p) {
        this.nodeMutate();
      }
      p -= 1;
    }
    p = this.mutationRates["enable"];
    while (p > 0) {
      if (Math.random() < p) {
        this.enableDisableMutation(true);
      }
      p -= 1;
    }
    p = this.mutationRates["disable"];
    while (p > 0) {
      if (Math.random() < p) {
        this.enableDisableMutation(false);
      }
      p -= 1;
    }
  }

  Genome.prototype.clone = function() {
    var genome = new Genome;
    for (var i = 0; i < this.genes.length; ++ i) {
      genome.genes.push(this.genes[i].clone());
    }
    genome.maxNeuron = this.maxNeuron;
    genome.mutationRates["connections"] = this.mutationRates["connections"];
    genome.mutationRates["link"] = this.mutationRates["link"];
    genome.mutationRates["bias"] = this.mutationRates["bias"];
    genome.mutationRates["node"] = this.mutationRates["node"];
    genome.mutationRates["enable"] = this.mutationRates["enable"];
    genome.mutationRates["disable"] = this.mutationRates["disable"];
    return genome;
  }

  Genome.disjoint = function(genome1, genome2) {
    var i1 = [];
    for (var i = 0; i < genome1.genes.length; ++ i) {
      var gene = genome1.genes[i];
      i1[gene.innovation] = true;
    }
    var i2 = [];
    for (var i = 0; i < genome2.genes.length; ++ i) {
      var gene = genome2.genes[i];
      i2[gene.innovation] = true;
    }
    var disjointGenes = 0;
    for (var i = 0; i < genome1.genes.length; ++ i) {
      var gene = genome1.genes[i];
      if (i2[gene.innovation] === undefined) {
        disjointGenes ++;
      }
    }
    for (var i = 0; i < genome2.genes.length; ++ i) {
      var gene = genome2.genes[i];
      if (i1[gene.innovation] === undefined) {
        disjointGenes ++;
      }
    }
    return disjointGenes / Math.max(genome1.genes.length, genome2.genes.length);
  }

  Genome.weights = function(genome1, genome2) {
    var i2 = [];
    for (var i = 0; i < genome2.genes.length; ++ i) {
      var gene = genome2.genes[i];
      i2[gene.innovation] = gene;
    }
    var sum = 0;
    var coincident = 0;
    for (var i = 0; i < genome1.genes.length; ++ i) {
      var gene = genome1.genes[i];
      if (i2[gene.innovation] !== undefined) {
        var gene2 = i2[gene.innovation];
        sum += Math.abs(gene.weight - gene2.weight);
        coincident ++;
      }
    }
    return sum / coincident;
  }

  Genome.sameSpecies = function(genome1, genome2) {
    var dd = DeltaDisjoint * Genome.disjoint(genome1, genome2);
    var dw = DeltaWeights * Genome.weights(genome1, genome2);
    return dd + dw < DeltaThreshold;
  }

  Genome.getBasicGenome = function() {
    var genome = new Genome();
    genome.maxNeuron = nInput;
    genome.mutate();
    return genome;
  }

  var Gene = function() {
    this.into = 0;
    this.out = 0;
    this.weight = 0;
    this.enabled = true;
    this.innovation = 0;
  }

  Gene.prototype.pack = function() {
    return JSON.stringify(this);
  }

  Gene.unpack = function(packed) {
    var storage = JSON.parse(packed);
    var g = new Gene;
    g.into = storage.into;
    g.out = storage.out;
    g.weight = storage.weight;
    g.enabled = storage.enabled;
    g.innovation = storage.innovation;
    return g;
  }

  Gene.prototype.clone = function() {
    var gene = new Gene();
    gene.into = this.into;
    gene.out = this.out;
    gene.weight = this.weight;
    gene.enabled = this.enabled;
    gene.innovation = this.innovation;
    return gene;
  }

  var Neuron = function() {
    // each input should be a gene
    this.inputs = [];
    this.value = 0;
  }

  Neuron.prototype.evaluate = function(neurons) {
    if (this.inputs.length > 0) {
      var sum = 0;
      for (var i = 0; i < this.inputs.length; ++ i) {
        var v = this.inputs[i].into;
        var w = this.inputs[i].weight;
        sum += neurons[v].value * w;
      }
      this.value = sigmoid(sum);
    }
  }

  // the first nInput neurons would be treated as inputs.
  // the order of neurons should be in topology order.
  var NeuronNetwork = function(genome) {
    this.neurons = [];
    for (var i = 0; i < nInput; ++ i) {
      this.neurons[i] = new Neuron();
    }
    for (var i = 0; i < nOutput; ++ i) {
      this.neurons[MaxNodes + i] = new Neuron();
    }
    genome.genes.sort(function (x, y) {
      return x.out - y.out;
    });

    for (var i = 0; i < genome.genes.length; ++ i) {
      var gene = genome.genes[i];
      if (gene.enabled === true) {
        if (this.neurons[gene.out] == undefined) {
          this.neurons[gene.out] = new Neuron();
        }

        var neuron = this.neurons[gene.out];
        neuron.inputs.push(gene);
        if (this.neurons[gene.into] === undefined) {
          this.neurons[gene.into] = new Neuron();
        }
      }
    }
    // genome.network = this;
  }

  NeuronNetwork.prototype.evaluate = function(inputs) {
    inputs.push(1); // add a constant one
    // assert (inputs.length() == nInput);

    for (var i = 0; i < inputs.length; ++ i) {
      this.neurons[i].value = inputs[i];
    }

    for (var idx in this.neurons) {
      this.neurons[idx].evaluate(this.neurons);
    }

    var outputs = {};

    for (var i = 0; i < nOutput; ++ i) {
      if (this.neurons[MaxNodes + i].value > 0) {
        outputs[inputKeys[i]] = true;
      } else {
        outputs[inputKeys[i]] = false;
      }
    }
    return outputs;
  }

  var NeatEngine = function(_game, _map, _brick, _bar) {
    this.pool = new Pool();
    this.pool.initialize();
    this.game = _game;
    this.map = _map;
    this.brick = _brick;
    this.bar = _bar;

    this.outputs = {};
  }

  NeatEngine.prototype.initializeRun = function() {
    $('#start_button').trigger("click");
    this.training = true;
    this.score = 0;
    this.pool.currentFrame = 0;
    this.timeout = TimeoutConstant;

    // clearJoypad()
    var sp = this.pool.species[this.pool.currentSpecies];
    var genome = sp.genomes[this.pool.currentGenome];
    genome.generateNetwork();
    this.evaluateCurrent();
  }

  NeatEngine.prototype.evaluateCurrent = function() {
    var sp = this.pool.species[this.pool.currentSpecies];
    var genome = sp.genomes[this.pool.currentGenome];
    
    var inputs = this.getInputs();
    var outputs = genome.network.evaluate(inputs);

    return outputs;
  }

  NeatEngine.prototype.nextGenome = function () {
    var sp = this.pool.species[this.pool.currentSpecies];
    this.pool.currentGenome ++;
    if (this.pool.currentGenome >= sp.genomes.length) {
      this.pool.currentGenome = 0;
      this.pool.currentSpecies ++;
      if (this.pool.currentSpecies >= this.pool.species.length) {
        this.pool.newGeneration();
        this.pool.currentSpecies = 0;
      }
    }
  }

  NeatEngine.prototype.fitnessAlreadyMeasured = function() {
    var sp = this.pool.species[this.pool.currentSpecies];
    var genome = sp.genomes[this.pool.currentGenome];
    return genome.fitness != 0;
  }

  NeatEngine.prototype.getInputs = function() {
    return readInputs(this.game, this.map, this.brick, this.bar);
  }

  NeatEngine.prototype.playTop = function() {
    var maxFitness = 0;
    var maxs, maxg;
    for (var i = 0; i < this.pool.species.length; ++ i) {
      var sp = this.pool.species[i];
      for (var j = 0; j < sp.genomes.length; ++ j) {
        var genome = sp.genomes[j];
        if (genome.fitness > maxFitness) {
          maxFitness = genome.fitness;
          maxs = i;
          maxg = j;
        }
      }
    }

    this.pool.currentSpecies = maxs;
    this.pool.currentGenome = maxg;
    this.pool.maxFitness = maxFitness;
    // console.log("Max Fitness: " + Math.floor(maxFitness));
    $('#maxFitness').html("Max Fitness: " + Math.floor(maxFitness));
    this.initializeRun();
    this.pool.currentFrame ++;
  }

  NeatEngine.prototype.sendOutputs = function() {
    for (var key in this.outputs) {
      if (this.outputs[key] == true) {
        if (this.nKeySent === undefined) {
          this.nKeySent = 0;
        }
        this.nKeySent ++;
        // console.log("send " + key);
        this.brick.handleInput(key);
      }
    }
  }

  NeatEngine.prototype.getScore = function() {
    var score = this.game.score * 100;
    
    for (var col = 0; col < NUM_COL; ++ col) {
      for (var row = 0; row < NUM_ROW; ++ row) {
        if (this.map.grid[col][row] == 3 || this.map.grid[col][row] == 4) {
          score += 50;
        }
        if (row < 5) {
          score += 10;
        }
      }
    }
    return score;
  }

  NeatEngine.prototype.takeOneStep = function() {
    var sp = this.pool.species[this.pool.currentSpecies];
    var genome = sp.genomes[this.pool.currentGenome];
    //if (this.pool.currentFrame % 5 == 0) {
      this.outputs = this.evaluateCurrent();
    //}
    this.sendOutputs();
    this.current_score = this.getScore();
    if (this.current_score > this.score) {
      this.score = this.current_score;
      this.timeout = TimeoutConstant;
    }

    this.timeout --;
    var timeoutBonus = this.pool.currentFrame / 4 + this.score / 4;
    if (this.timeout + timeoutBonus <= 0 || !this.game.status) {
      var fitness = this.score - this.pool.currentFrame / 2;
      if (fitness == 0) fitness = -1;
      genome.fitness = fitness;
      if (fitness > this.pool.maxFitness) {
        this.pool.maxFitness = fitness;
        //console.log("Max Fitness: " + Math.floor(fitness));
        $('#maxFitness').html("Max Fitness: " + Math.floor(fitness));
      }
      //console.log("Gen " + this.pool.generation + " sp " + this.pool.currentSpecies
      //    + " genome " + this.pool.currentGenome + " fitness " + fitness);
      $('#neat_msg').html("Gen " + this.pool.generation + " sp " +
          this.pool.currentSpecies + " genome " + this.pool.currentGenome +
          " fitness " + fitness);
      this.pool.currentSpecies = 0;
      this.pool.currentGenome = 0;
      while (this.fitnessAlreadyMeasured()) {
        this.nextGenome();
      }
      this.initializeRun();
    }

    var measured = 0;
    var total = 0;

    for (var i = 0; i < this.pool.species.length; ++ i) {
      var sp = this.pool.species[i];
      for (var j = 0; j < sp.genomes.length; ++ j) {
        total ++;
        if (sp.genomes[j].fitness != 0) {
          measured ++;
        }
      }
    }
    //console.log("Gen " + this.pool.generation + " measured: " + measured + "/" + total);
    $('#neat_msg').html("Gen " + this.pool.generation + " measured: " + measured + "/" + total);
    this.pool.currentFrame ++;
  }

  var neat;
  var tid;

  function startAI(interval, game, map, brick, bar) {
    console.log("starting AI");
    // $('#start_button').trigger("click");
    if (neat === undefined) {
      neat = new NeatEngine(game, map, brick, bar);
    }
    neat.initializeRun();
    // tid = setInterval(function () { neat.takeOneStep(); }, interval);
    // window.URL.revokeObjectURL
    // window.URL.createObjectURL
  }

  function stopAI() {
    if (tid !== undefined) clearInterval(tid);
    tid = undefined;
    neat.training = false;
    // neat = undefined;
  }

  $(document).ready(function() {
    $('#neat_run').click(function() {
      startAI(100, game, map, brick, bar);
    });
    $('#neat_stop').click(function() {
      stopAI();
    });
  });
/* })(); */
