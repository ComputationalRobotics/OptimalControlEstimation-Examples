import torch
import torch.nn as nn
import torch.nn.init as init

class CostToGo(nn.Module):
    def __init__(self, input_size=2, hidden_size=3072, output_size=1024):
        super(CostToGo, self).__init__()
        self.fc1 = nn.Linear(input_size, hidden_size, bias=True)
        self.fc2 = nn.Linear(hidden_size, hidden_size, bias=True)
        self.fc3 = nn.Linear(hidden_size, output_size, bias=True)
        self.activation = nn.ReLU()

    def forward(self, x):
        x = self.activation(self.fc1(x))
        x = self.activation(self.fc2(x))
        x = self.fc3(x)
        return x