#!/usr/bin/env python3
"""
Script for a simple neural netowrk to predict phenotype from hog pathway P(Aff) scores.
Returns column names (i.e. gene IDs), paff scores (features), S-scores (labels)
"""
import pandas as pd
import scipy
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
import tensorflow as tf
from tensorflow.keras import layers

MODEL_DIR = '/Users/ally/Projects/predictor/models'
HOG_PAFF_PATH = '/Users/ally/Projects/predictor/data/hog_paff.tsv'
HOG_GROWTH_PATH = '/Users/ally/Projects/predictor/data/liti_yeast_growth.tsv'

BATCH_SIZE = 25

def import_hog_data(paff_path, growth_path):
    """
    Load Phenotype and HOG pathway P(Aff) scores into a combined pd.DataFrame
    """
    growth = pd.read_csv(growth_path, sep='\t', usecols=['strain', 'ypdnacl15m']).dropna(0)
    paff = pd.read_csv(paff_path, sep='\t')
    #paff = paff[paff['strain'].isin(growth['strain'])]

    paff = paff.set_index('strain')
    growth = growth.set_index('strain')

    return pd.concat([paff, growth], axis=1, join='inner')

def paff_dataset(features, labels, batch_size, shuffle_buffer=0, repeat=False):
    """
    Create tf.Dataset from HOG pathway and phenotype data
    """
    data = (tf.cast(features, tf.float32), tf.cast(labels, tf.float32))
    dataset = tf.data.Dataset.from_tensor_slices(data)

    if shuffle_buffer:
        dataset = dataset.shuffle(shuffle_buffer)

    if repeat:
        dataset = dataset.repeat()

    dataset = dataset.batch(batch_size)
    return dataset

def create_model(layer_sizes, drop_rate):
    """
    Create a HOG predictor neural network
    """
    model = tf.keras.Sequential()
    for i in layer_sizes:
        model.add(layers.Dense(i))
        model.add(layers.Dropout(drop_rate))

    model.add(layers.Dense(1))

    return model

def main():
    """
    Main function
    """
    # Load data
    hog_data = import_hog_data(HOG_PAFF_PATH, HOG_GROWTH_PATH)

    feature_columns = hog_data.columns.values[:-1]
    label_column = hog_data.columns.values[-1]

    train, test = train_test_split(hog_data, test_size=0.1,
                                   shuffle=True, random_state=723534)

    training_dataset = paff_dataset(train[feature_columns], train[label_column],
                                    batch_size=BATCH_SIZE, shuffle_buffer=75, repeat=True)

    test_dataset = paff_dataset(test[feature_columns], test[label_column],
                                batch_size=len(test), repeat=True)

    model = create_model([100, 100, 50], 0.1)

    model.compile(loss='mean_absolute_error',
                  optimizer=tf.keras.optimizers.SGD(lr=0.1,
                                                    momentum=0.01,
                                                    clipvalue=0.1),
                  metrics=['mean_absolute_error'])

    ktb = tf.keras.callbacks.TensorBoard(log_dir=f'{MODEL_DIR}/hog_predictor',
                                         histogram_freq=5, write_graph=True)

    name = 'chkpt_{epoch:02d}'
    kcp = tf.keras.callbacks.ModelCheckpoint(filepath=f"{MODEL_DIR}/hog_predictor/{name}",
                                             period=25)

    # Determine number of steps (e.g. number of minibatches per epoch)
    train_steps = len(train) // BATCH_SIZE + 1
    val_steps = 1 # set batch to include all

    model.fit(training_dataset,
              epochs=50,
              steps_per_epoch=train_steps,
              validation_data=test_dataset,
              callbacks=[ktb, kcp],
              validation_steps=val_steps)

    model.save(f"{MODEL_DIR}/hog_predictor/final_model.h5")

    hog_data['nn_pred'] = model.predict(hog_data[feature_columns])
    hog_data['test'] = 0
    hog_data.loc[hog_data.index.intersection(test.index), 'test'] = 1

    # Compare to linear model
    linear_model = LinearRegression()
    linear_model.fit(hog_data[feature_columns], hog_data[label_column])
    hog_data['lm_pred'] = linear_model.predict(hog_data[feature_columns])

    cols = ['r' if i else 'k' for i in hog_data['test']]
    fig, axes = plt.subplots(2, 2, figsize=(9, 9))
    fig.suptitle('Predicting yeast growth in 1.5mM NaCl from HOG P(aff) scores')

    axes[0, 0].set_title('NN Predictor')
    axes[0, 0].set_xlabel('True Growth')
    axes[0, 0].set_ylabel('Predicted Growth')
    axes[0, 0].scatter(hog_data[label_column], hog_data['nn_pred'], c=cols)

    axes[0, 1].set_title('Linear Model')
    axes[0, 1].set_xlabel('True Growth')
    axes[0, 1].set_ylabel('Predicted Growth')
    axes[0, 1].scatter(hog_data[label_column], hog_data['lm_pred'], c=cols)

    axes[1, 0].set_title('Comparison 1')
    axes[1, 0].set_xlabel('Linear Model')
    axes[1, 0].set_ylabel('Neural Network')
    axes[1, 0].scatter(hog_data['lm_pred'], hog_data['nn_pred'], c=cols)

    axes[1, 1].set_title('Comparison 2')
    axes[1, 1].set_xlabel('True Growth')
    axes[1, 1].set_ylabel('Predicted Growth')
    axes[1, 1].scatter(hog_data[label_column], hog_data['nn_pred'], c='b')
    axes[1, 1].scatter(hog_data[label_column], hog_data['lm_pred'], c='g')

    patches = [mpatches.Patch(color='r', label='Test'),
               mpatches.Patch(color='k', label='Train'),
               mpatches.Patch(color='b', label='NN'),
               mpatches.Patch(color='g', label='LM')]
    plt.legend(handles=patches)
    fig.savefig('/Users/ally/Projects/predictor/figures/hog_predictor.pdf')

    # Correlations
    print('All data:')
    lm_cor, lm_p = scipy.stats.pearsonr(hog_data[label_column], hog_data['lm_pred'])
    nn_cor, nn_p = scipy.stats.pearsonr(hog_data[label_column], hog_data['nn_pred'])
    print(f'\tLinear model: r = {lm_cor}, p = {lm_p}')
    print(f'\tNeural Network: r = {nn_cor}, p = {nn_p}')

    print('Training set:')
    lm_cor, lm_p = scipy.stats.pearsonr(hog_data.loc[hog_data.test == 0, label_column],
                                        hog_data.loc[hog_data.test == 0, 'lm_pred'])

    nn_cor, nn_p = scipy.stats.pearsonr(hog_data.loc[hog_data.test == 0, label_column],
                                        hog_data.loc[hog_data.test == 0, 'nn_pred'])
    print(f'\tLinear model: r = {lm_cor}, p = {lm_p}')
    print(f'\tNeural Network: r = {nn_cor}, p = {nn_p}')

    print('Test set:')
    lm_cor, lm_p = scipy.stats.pearsonr(hog_data.loc[hog_data.test == 1, label_column],
                                        hog_data.loc[hog_data.test == 1, 'lm_pred'])

    nn_cor, nn_p = scipy.stats.pearsonr(hog_data.loc[hog_data.test == 1, label_column],
                                        hog_data.loc[hog_data.test == 1, 'nn_pred'])
    print(f'\tLinear model: r = {lm_cor}, p = {lm_p}')
    print(f'\tNeural Network: r = {nn_cor}, p = {nn_p}')

if __name__ == "__main__":
    main()
