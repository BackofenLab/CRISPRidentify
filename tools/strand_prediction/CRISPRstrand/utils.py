import numpy as np
#import matplotlib.pyplot as plt
from keras.models import load_model as lm

# =============================================================================
# 
# =============================================================================
def load_model(model_path):
    classifier = lm(model_path)
    classifier.summary()
    sequence_length = classifier.layers[1].input_shape[2]

    return classifier, sequence_length

def plot_loss_acc(histories):
    
    training_accuracies = []
    validation_accuracies = []
    training_losses = []
    validation_losses = []
    
    for history in histories:
        training_accuracies.append(history.history['acc'])
        validation_accuracies.append(history.history['val_acc'])
        training_losses.append(history.history['loss'])
        validation_losses.append(history.history['val_loss'])
    
    training_accuracies = np.stack(training_accuracies)
    validation_accuracies = np.stack(validation_accuracies)
    training_losses = np.stack(training_losses)
    validation_losses = np.stack(validation_losses)

    mean_training_accuracies = np.mean(training_accuracies, axis = 0)
    mean_validation_accuracies = np.mean(validation_accuracies, axis = 0)
    mean_training_losses = np.mean(training_losses, axis = 0)
    mean_validation_losses = np.mean(validation_losses, axis = 0)
    
    std_training_accuracies = np.sqrt(np.var(training_accuracies, axis = 0))
    std_validation_accuracies = np.sqrt(np.var(validation_accuracies, axis = 0))
    std_training_losses = np.sqrt(np.var(training_losses, axis = 0))
    std_validation_losses = np.sqrt(np.var(validation_losses, axis = 0))
    
    epochs = training_accuracies.shape[1]
    
    x = np.atleast_2d(np.linspace(1, epochs, epochs)).T
    plt.plot(x, mean_training_accuracies, label='train_acc')
    plt.plot(x, mean_validation_accuracies, label='val_acc')
    plt.fill(np.concatenate([x, x[::-1]]), np.concatenate([mean_training_accuracies-std_training_accuracies, (mean_training_accuracies+std_training_accuracies)[::-1]]), alpha=.3, label='train_acc_uncertainty')
    plt.fill(np.concatenate([x, x[::-1]]), np.concatenate([mean_validation_accuracies-std_validation_accuracies, (mean_validation_accuracies+std_validation_accuracies)[::-1]]), alpha=.3, label='val_acc_uncertainty')
    plt.legend(loc='lower right')
    plt.show()
    
    plt.plot(x, mean_training_losses, label='train_loss')
    plt.plot(x, mean_validation_losses, label='val_loss')
    plt.fill(np.concatenate([x, x[::-1]]), np.concatenate([mean_training_losses-std_training_losses, (mean_training_losses+std_training_losses)[::-1]]), alpha=.3, label='train_loss_uncertainty')
    plt.fill(np.concatenate([x, x[::-1]]), np.concatenate([mean_validation_losses-std_validation_losses, (mean_validation_losses+std_validation_losses)[::-1]]), alpha=.3, label='val_loss_uncertainty')
    plt.legend(loc='upper right')
    plt.show()
