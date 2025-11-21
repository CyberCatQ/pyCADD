import numpy as np
from pandas import DataFrame
import torch
from torch import nn


class MLP(nn.Module):
    """Multi-Layer Perceptron implemented with PyTorch.

    A fully connected neural network for molecular property prediction
    and classification tasks in CADD workflows.
    """

    def __init__(
        self, input_dim: int, hidden_dim: list, output_dim: int, device: torch.device = None
    ) -> None:
        """Initialize the MLP model.

        Args:
            input_dim (int): Dimension of input features.
            hidden_dim (list): List of dimensions for hidden layers.
            output_dim (int): Dimension of output layer.
            device (torch.device, optional): Device to run the model on (CPU or GPU).
        """
        super(MLP, self).__init__()
        self.input_dim = input_dim
        self.hidden_dim = hidden_dim
        self.output_dim = output_dim
        self.device = (
            device
            if device is not None
            else torch.device("cuda" if torch.cuda.is_available() else "cpu")
        )

        # Number of hidden layers
        self.layer_num = len(hidden_dim)

        # Define model architecture
        model = nn.Sequential()
        model.add_module("input_layer", nn.Linear(input_dim, hidden_dim[0]))
        model.add_module("relu_1", nn.ReLU())

        for layer_idx in range(self.layer_num - 1):
            model.add_module(
                f"hidden_layer_{layer_idx + 1}",
                nn.Linear(hidden_dim[layer_idx], hidden_dim[layer_idx + 1]),
            )
            model.add_module(f"relu_{layer_idx + 2}", nn.ReLU())

        model.add_module("output_layer", nn.Linear(hidden_dim[-1], output_dim))

        self.model = model.to(device)

    def forward(self, x):
        return self.model(x)

    def get_params(self) -> dict:
        """Get model parameters.

        Returns:
            Dictionary containing model configuration parameters.
        """
        return {
            "device": self.device,
            "input_dim": self.input_dim,
            "output_dim": self.output_dim,
            "hidden_dim": self.hidden_dim,
        }

    def initialize(self, method: str = "normal") -> None:
        """Initialize model parameters.

        Args:
            method (str): Initialization method. Options:
                - 'normal': Normal distribution initialization
                - 'xavier_uniform': Xavier uniform initialization
                - 'xavier_normal': Xavier normal distribution initialization
                - 'kaiming_uniform': Kaiming uniform initialization
                - 'kaiming_normal': Kaiming normal distribution initialization

        Raises:
            ValueError: If initialization method is unknown.
        """
        if method == "xavier_uniform":
            init_func = nn.init.xavier_uniform_
        elif method == "xavier_normal":
            init_func = nn.init.xavier_normal_
        elif method == "normal":
            init_func = nn.init.normal_
        elif method == "kaiming_uniform":
            init_func = nn.init.kaiming_uniform_
        elif method == "kaiming_normal":
            init_func = nn.init.kaiming_normal_
        else:
            raise ValueError("Unknown initialization method")
        for layer in self.model:
            if isinstance(layer, nn.Linear):
                init_func(layer.weight)

    def predict(self, X: DataFrame) -> np.ndarray:
        """Make predictions on input data.

        Args:
            X (DataFrame): Input feature data.

        Returns:
            Numpy array containing predicted class labels.
        """
        self.eval()
        X = torch.tensor(X.values, dtype=torch.float)
        with torch.no_grad():
            X = X.to(self.device)
            y = self.forward(X)
            y = y.cpu().argmax(dim=1).numpy()
        return y

    def predict_proba(self, X: DataFrame) -> np.ndarray:
        """Make probability predictions on input data.

        Args:
            X (DataFrame): Input feature data.

        Returns:
            Numpy array containing predicted class probabilities for each class.
        """
        self.eval()
        X = torch.tensor(X.values, dtype=torch.float)
        with torch.no_grad():
            X = X.to(self.device)
            y = self.forward(X)
            y = y.cpu().softmax(dim=1).numpy()
        return y
