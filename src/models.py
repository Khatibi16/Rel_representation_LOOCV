import torch
import torch.nn as nn

class DeepSurv(nn.Module):
    def __init__(self, input_dim,
                 enc_dims,
                 surv_dims,
                 proj_dims,
                 dropout):
        super(DeepSurv, self).__init__()
        self.n_layers_surv = len(surv_dims)
        n_layers_enc =len(enc_dims)

        #init encoder
        self.encoder = nn.ModuleList([nn.Sequential(
            nn.Linear(input_dim, enc_dims[0]),
            nn.ReLU(),
            nn.BatchNorm1d(enc_dims[0]),
            nn.Dropout(dropout)
        )])

        for layer in range(1, n_layers_enc):
            self.encoder.append(nn.Sequential(
                nn.Linear(enc_dims[layer - 1], enc_dims[layer]),
                nn.ReLU(),
                nn.BatchNorm1d(enc_dims[layer]), #todo: check what happens during transfer
                nn.Dropout(dropout)
            ))

        #init survival head
        if self.n_layers_surv > 0:
            #init surv module
            self.surv_module = nn.ModuleList([nn.Sequential(
                nn.Linear(enc_dims[-1], surv_dims[0]),
                nn.ReLU(),
                nn.BatchNorm1d(surv_dims[0]),
                nn.Dropout(dropout)
            )])

            for layer in range(1, self.n_layers_surv):
                self.surv_module.append(nn.Sequential(
                    nn.Linear(surv_dims[layer - 1], surv_dims[layer]),
                    nn.ReLU(),
                    nn.BatchNorm1d(surv_dims[layer]),
                    nn.Dropout(dropout)
                ))

            self.output = nn.Linear(surv_dims[-1], 1)
        else:
            self.output = nn.Linear(enc_dims[-1], 1)

        #init projection head
        self.proj_head=nn.Linear(enc_dims[-1], proj_dims)

    def forward(self, x, get_projection=False, get_emb = False):
        for layer in self.encoder:
            x = layer(x)

        x_proj = self.proj_head(x) #todo potentially remove

        if self.n_layers_surv > 0:
            for layer in self.surv_module:
                x= layer(x)

        if get_projection:
            return x_proj
        elif get_emb:
            return x_proj.detach()
        else:
            return self.output(x)

    def freeze_encoder_layers(self, n_layers):
        """
        Freezes the first `n_layers` of the encoder.
        """
        for i, layer in enumerate(self.encoder):
            if i < n_layers:
                for param in layer.parameters():
                    param.requires_grad = False
