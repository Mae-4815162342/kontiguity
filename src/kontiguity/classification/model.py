"""Built-in classifier, stored as plain JSON weights (no pickle, no scikit-learn needed)."""
import json
import os

import numpy as np

BUILTIN = os.path.join(os.path.dirname(__file__), "models", "hitchhiker_mlp.json")


class MLPEnsemble:
    """Mean probability of several (StandardScaler -> dense ReLU layers -> logistic) networks."""

    def __init__(self, path="builtin"):
        path = BUILTIN if path in (None, "", "builtin") else path
        with open(path) as fh:
            self.spec = json.load(fh)
        self.features = self.spec["features"]
        self.resolution = int(self.spec.get("resolution", 5000))

    @staticmethod
    def _forward(m, X):
        h = (X - np.asarray(m["mean"])) / np.asarray(m["scale"])
        W, B = m["weights"], m["biases"]
        for i, (w, b) in enumerate(zip(W, B)):
            h = h @ np.asarray(w) + np.asarray(b)
            if i < len(W) - 1:
                h = np.maximum(h, 0) if m.get("activation", "relu") == "relu" else np.tanh(h)
        return 1 / (1 + np.exp(-h[:, 0]))

    def predict_proba(self, X):
        X = np.atleast_2d(np.asarray(X, float))
        return np.mean([self._forward(m, X) for m in self.spec["members"]], 0)
