
class SIR_node():

    def __init__(self):
        self.set_susceptible()

    def is_infected(self):
        return self.infected

    def is_recovered(self):
        return self.recovered

    def is_susceptible(self):
        return not self.infected and not self.recovered

    def set_infected(self):
        self.infected = True
        self.recovered = False

    def set_recovered(self):
        self.infected = False
        self.recovered = True

    def set_susceptible(self):
        self.infected = False
        self.recovered = False


class SEIR_node():

    def __init__(self):
        self.set_susceptible()

    def is_infected(self):
        return self.infected

    def is_recovered(self):
        return self.recovered

    def is_exposed(self):
        return self.exposed

    def is_susceptible(self):
        return not self.infected and not self.recovered and not self.exposed

    def set_infected(self):
        self.infected = True
        self.recovered = False
        self.exposed = False

    def set_exposed(self):
        self.infected = False
        self.recovered = False
        self.exposed = True

    def set_recovered(self):
        self.infected = False
        self.recovered = True
        self.exposed = False

    def set_susceptible(self):
        self.infected = False
        self.recovered = False
        self.exposed = False
