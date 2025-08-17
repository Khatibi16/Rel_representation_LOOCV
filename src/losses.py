
import numpy as np
import torch

def risk_set_matrix(dim):
    # create indicator matrix
    risk_sets = np.zeros((dim, dim))
    il = np.tril_indices(len(risk_sets))
    risk_sets[il] = 1
    risk_sets = torch.from_numpy(risk_sets)

    return risk_sets

def nlplloss_vectorized( event, time, hazard, device, per_sample=False):
    # create a matrix of time, event, and hazard
    concat = torch.cat((time.reshape(-1, 1), event.reshape(-1, 1), hazard.reshape(-1, 1)), dim=1)  # changed axis to dim

    # sort matrix according to survival time
    concat = concat[concat[:, 0].sort()[1]]

    # split into event and hazard
    time = concat[:, 1].reshape(-1, 1)
    event = concat[:, 1].reshape(-1, 1)
    hazard = concat[:, 2].reshape(-1, 1)
    del concat

    n_events = torch.sum(event)

    # getting risk sets
    risk_sets = risk_set_matrix(hazard.shape[0]).to(device)

    # calculate risk set hazards
    risk_set_haz = torch.mul(torch.exp(hazard), risk_sets)

    # sum the hazards over the risk sets
    summed_risk_set_haz = torch.log(torch.sum(risk_set_haz, dim=0) + 1e-4)
    summed_risk_set_haz = torch.reshape(summed_risk_set_haz, (-1, 1))

    # subtract risk set hazards from individual hazards
    subtracted_haz = torch.reshape(torch.sub(hazard, summed_risk_set_haz), (-1, 1))

    if per_sample:
        # count only those with events for final loss
        loss = torch.mul(subtracted_haz, event).float()

    else:
        # count only those with events for final loss
        loss = torch.sum(torch.mul(subtracted_haz, event))

        # normalize by mini batch
        loss = torch.div(loss, n_events + 1)

    return -loss