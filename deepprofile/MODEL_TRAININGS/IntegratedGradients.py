import numpy as np
import tensorflow as tf

class IntegratedGradients:
    def __init__(self, model):
        self.model = model

    def interpolate_inputs(self, baseline, input, steps):
        alphas = np.linspace(0, 1, steps)
        delta = input - baseline
        interpolated = np.array([baseline + alpha * delta for alpha in alphas])
        return interpolated

    def compute_gradients(self, inputs, target_neuron):
        with tf.GradientTape() as tape:
            tape.watch(inputs)
            outputs = self.model(inputs)
            if len(outputs.shape) == 2:
                out = outputs[:, target_neuron]
            else:
                out = outputs[target_neuron]
        grads = tape.gradient(out, inputs)
        return grads.numpy()

    def explain(self, input, target_neuron=0, baseline=None, steps=50):
        if baseline is None:
            baseline = np.zeros_like(input)
        interpolated_inputs = self.interpolate_inputs(baseline, input, steps)
        grads = []
        for x in interpolated_inputs:
            x_tf = tf.convert_to_tensor(np.expand_dims(x, axis=0), dtype=tf.float32)
            grad = self.compute_gradients(x_tf, target_neuron)
            grads.append(grad[0])
        avg_grads = np.mean(grads, axis=0)
        integrated_gradients = (input - baseline) * avg_grads
        return integrated_gradients
