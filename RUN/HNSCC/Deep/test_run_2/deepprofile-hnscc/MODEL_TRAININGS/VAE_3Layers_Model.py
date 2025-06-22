import os
import numpy as np
import pandas as pd
from sklearn.metrics import mean_squared_error, r2_score
import tensorflow as tf
from tensorflow.keras.layers import (Input, Dense, Activation, BatchNormalization, Layer)
from tensorflow.keras.models import Model
from tensorflow.keras import optimizers
from tensorflow.keras.callbacks import Callback
import sys

# --- Args ---
cancer_type = sys.argv[1]
intermediate1_dim = int(sys.argv[2])
intermediate2_dim = int(sys.argv[3])
latent_dim = int(sys.argv[4])
fold = int(sys.argv[5])

input_folder = './ALL_CANCER_FILES/' + cancer_type + '/'
output_folder = './ALL_CANCER_FILES/' + cancer_type + '/VAE_FILES/'
os.makedirs(output_folder, exist_ok=True)

input_filename = input_folder + cancer_type + '_DATA_TOP2_JOINED_PCA_500L.tsv'
output_filename = cancer_type + '_DATA_TOP2_JOINED_encoded_'

input_df = pd.read_table(input_filename, index_col=0)
print("INPUT FILE", input_df.shape)
print(input_df.head(5))

original_dim = input_df.shape[1]
batch_size = 50
epochs = 50
learning_rate = 0.0005
beta = tf.Variable(1.0)
kappa = 0

np.random.seed(123456 * fold)
tf.random.set_seed(123456 * fold)

# ---- Modern Sampling Layer with KL Loss ----
class Sampling(Layer):
    def __init__(self, beta, **kwargs):
        super().__init__(**kwargs)
        self.beta = beta

    def call(self, inputs):
        z_mean, z_log_var = inputs
        epsilon = tf.random.normal(shape=tf.shape(z_mean))
        z = z_mean + tf.exp(0.5 * z_log_var) * epsilon
        kl_loss = -0.5 * tf.reduce_sum(1 + z_log_var - tf.square(z_mean) - tf.exp(z_log_var), axis=1)
        self.add_loss(self.beta * tf.reduce_mean(kl_loss))
        return z

# ---- Encoder ----
x = Input(shape=(original_dim,))
net = Dense(intermediate1_dim, kernel_initializer='glorot_uniform')(x)
net = BatchNormalization()(net)
net = Activation('relu')(net)
net = Dense(intermediate2_dim, kernel_initializer='glorot_uniform')(net)
net = BatchNormalization()(net)
net = Activation('relu')(net)
z_mean = Dense(latent_dim, kernel_initializer='glorot_uniform', name="z_mean")(net)
z_log_var = Dense(latent_dim, kernel_initializer='glorot_uniform', name="z_log_var")(net)
z = Sampling(beta)([z_mean, z_log_var])

encoder = Model(x, z_mean, name="encoder")
encoder.summary()

# ---- Decoder ----
latent_inputs = Input(shape=(latent_dim,))
h_decoded = Dense(intermediate2_dim, activation='relu', kernel_initializer='glorot_uniform')(latent_inputs)
h_decoded = Dense(intermediate1_dim, activation='relu', kernel_initializer='glorot_uniform')(h_decoded)
x_decoded_mean = Dense(original_dim, kernel_initializer='glorot_uniform')(h_decoded)
decoder = Model(latent_inputs, x_decoded_mean, name="decoder")
decoder.summary()

# ---- VAE Model ----
outputs = decoder(z)
vae = Model(x, outputs, name="vae")

optimizer = optimizers.Adam(learning_rate=learning_rate)
vae.compile(optimizer=optimizer, loss="mse")
vae.summary()

# ---- Beta Warm-up Callback ----
class WarmUpCallback(Callback):
    def __init__(self, beta, kappa):
        super().__init__()
        self.beta = beta
        self.kappa = kappa

    def on_epoch_end(self, epoch, logs=None):
        if self.beta.numpy() <= 1:
            self.beta.assign(self.beta.numpy() + self.kappa)

# ---- Train ----
history = vae.fit(
    np.array(input_df), np.array(input_df),
    shuffle=True,
    epochs=epochs,
    batch_size=batch_size,
    verbose=2,
    callbacks=[WarmUpCallback(beta, kappa)]
)

# ---- Encode, reconstruct, and evaluate ----
training_encoded = encoder.predict(np.array(input_df), batch_size=batch_size)
training_encoded_df = pd.DataFrame(training_encoded, index=input_df.index)

training_reconstructed = decoder.predict(training_encoded, batch_size=batch_size)
training_reconstructed_df = pd.DataFrame(training_reconstructed, index=input_df.index, columns=input_df.columns)

recons_error = mean_squared_error(np.array(input_df), np.array(training_reconstructed_df))
print("TRAINING RECONSTRUCTION ERROR: " + str(recons_error))

training_encoded_df.to_csv(
    output_folder + output_filename + str(latent_dim) + "L_TRAINING_fold" + str(fold) + ".tsv",
    sep='\t'
)

encoder_json = encoder.to_json()
with open(output_folder + f"VAE_{cancer_type}_encoder_{latent_dim}L_{fold}.json", "w") as json_file:
    json_file.write(encoder_json)
encoder.save_weights(output_folder + f"VAE_{cancer_type}_encoder_{latent_dim}L_{fold}.weights.h5")
print("Saved encoder model to disk")

decoder_json = decoder.to_json()
with open(output_folder + f"VAE_{cancer_type}_decoder_{latent_dim}L_{fold}.json", "w") as json_file:
    json_file.write(decoder_json)
decoder.save_weights(output_folder + f"VAE_{cancer_type}_decoder_{latent_dim}L_{fold}.weights.h5")
print("Saved decoder model to disk")

training_r2_vals = np.zeros(input_df.shape[0])
for i in range(input_df.shape[0]):
    training_r2 = r2_score(input_df.values[i, :], training_reconstructed_df.values[i, :])
    training_r2_vals[i] = training_r2
print("TRAINING R2 " + str(np.mean(training_r2_vals)))
