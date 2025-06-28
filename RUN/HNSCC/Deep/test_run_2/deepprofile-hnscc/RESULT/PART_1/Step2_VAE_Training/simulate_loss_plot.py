# simulate_loss_plot.py
import matplotlib.pyplot as plt
import numpy as np

# Simulate
epochs = list(range(1, 101))
train_loss = np.exp(-0.05 * np.array(epochs)) + 0.01 * np.random.rand(100)
val_loss = np.exp(-0.045 * np.array(epochs)) + 0.02 * np.random.rand(100)

# Plot
plt.figure(figsize=(8, 5))
plt.plot(epochs, train_loss, label="Training Loss")
plt.plot(epochs, val_loss, label="Validation Loss")
plt.xlabel("Epoch")
plt.ylabel("Loss")
plt.title("Simulated VAE Training Loss")
plt.legend()
plt.tight_layout()
plt.savefig("vae_loss_curve.png")
plt.savefig("vae_loss_curve.pdf")
plt.show()
