import torch as th
import numpy as np
from torch.utils.data import DataLoader, TensorDataset
import torch.optim as optim
import control as ctl
from tqdm import tqdm
from typing import Dict, List, Tuple
from torch.func import vmap, jacrev, functional_call
import utils
from CostToGoFunction import CostToGo
from LinearDynamics import LinearSystem

# Device and Seed Initialization
device = th.device("cuda" if th.cuda.is_available() else "cpu")
seed = th.randint(0, 10000000, (1,))
th.manual_seed(seed)
if th.cuda.is_available():
    th.cuda.manual_seed_all(seed)
np.random.seed(seed)
print(f"seed = {seed}")

########################################################################
############################## PARAMETERS ##############################
########################################################################

batch_size = 32768
epoch = 10000
U_MAX = float('inf')
x0max = 8
x1max = 8
batch_size_train = 2048
pretrain = True

model = CostToGo(input_size=2).to(device)
if pretrain:
    model.load_state_dict(th.load('./CTOptimalControlExperiment/Linear/model/initial_model.pth'))
optimizer = optim.Adam(model.parameters(), lr=0.00005)
scheduler = StepLR(optimizer, step_size=50, gamma=0.99)
system = LinearSystem()

def functionalized_model(
    model: CostToGo,
    params: Dict,
    x: th.Tensor,
):  
    res = functional_call(model, params, x)
    return res @ th.t(res)

def get_fx_gx_torch(
    x: th.Tensor,
    SystemModel: system,
) -> Tuple[th.Tensor, th.Tensor]:
    
    batch_size, x_dim = x.shape
    fx, gx = system.dynamics_torch(x, device)

    fx = fx.unsqueeze(2)
    gx = gx.unsqueeze(2)
    return fx, gx

def LQR_controller(xt):
    A = th.tensor([[0, 1],
              [0, 0]], dtype=th.float32, device=device)

    B = th.tensor([[0],
                [1]], dtype=th.float32, device=device)

    Q = th.tensor([[1, 0],
                [0, 1]], dtype=th.float32, device=device)

    R = th.tensor([[1]], dtype=th.float32, device=device)

    K, S, E = ctl.lqr(A.cpu(), B.cpu(), Q.cpu(), R.cpu())
    th_K = th.tensor(K, dtype=th.float32, device=device)
    u = -th_K @ xt
        
    return u
    

if __name__ == "__main__":
    x = th.tensor(utils.prepare_data(sample_number=batch_size, x0max=x0max, x1max=x1max, vector_length=2), dtype=th.float32, requires_grad=True).to(device)
    x_eval = th.tensor(utils.prepare_data(sample_number=16 * batch_size_train, x0max=x0max, x1max=x1max, vector_length=2), dtype=th.float32, requires_grad=True).to(device)
    
    # DataLoader Initialization
    train_dataset1 = TensorDataset(x)
    train_loader1 = DataLoader(train_dataset1, batch_size=batch_size_train, shuffle=True)
    eval_dataset1 = TensorDataset(x_eval)
    eval_loader1 = DataLoader(eval_dataset1, batch_size=batch_size_train, shuffle=True)

    single_dJdx_func = jacrev(functionalized_model, argnums=2)
    batch_dJdx_func = vmap(single_dJdx_func, in_dims=(None, None, 0))
    batch_th_dot = vmap(th.dot)
    batch_LQR_controller = vmap(LQR_controller)

    # loss 
    loss_values = []
    best_val_loss = float('inf')

    # Training
    for i in range(epoch):
        for (batch_x,)in tqdm(train_loader1):
            batch_dJdx = batch_dJdx_func(model, dict(model.named_parameters()), batch_x)
            batch_lqr_u = th.squeeze(batch_LQR_controller(batch_x))
            fx, gx = get_fx_gx_torch(batch_x, system)
            batch_dJdx = batch_dJdx.unsqueeze(1)

            LfJ = th.squeeze(batch_dJdx @ fx) # [batch_size, ]
            LgJ = th.squeeze(batch_dJdx @ gx) # [batch_size, ]
            batch_x_square = batch_th_dot(batch_x, batch_x) # [batch_size, ]
            
            batch_h1 = (batch_x_square + LfJ - 0.25 * (LgJ**2))    # [batch_size, ]
            
            sign_loss = 0.05*th.sum(th.relu(LgJ * batch_lqr_u))
            wrong_number = 0.5*th.sum(th.abs(th.sign(LgJ) + th.sign(batch_lqr_u)))
            assert th.all(th.relu(LgJ * batch_lqr_u).ge(0))
            
            x0 = th.tensor([0.0,0.0], dtype=th.float32, device=device)
            residual_loss = (th.sum(th.abs(batch_h1)))/batch_size_train + th.dot(model(x0), model(x0))
            loss = sign_loss + residual_loss
            
            loss_values.append(loss.item())
            
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            scheduler.step()
            
        model.eval()
        with th.no_grad():
            loss_eval = 0
            for (batch_x_eval,) in tqdm(eval_loader1):
                batch_dJdx_eval = batch_dJdx_func(model, dict(model.named_parameters()), batch_x_eval)
                fx_eval, gx_eval = get_fx_gx_torch(batch_x_eval, system)
                batch_lqr_u_eval = th.squeeze(batch_LQR_controller(batch_x_eval))
                
                batch_dJdx_eval = batch_dJdx_eval.unsqueeze(1)

                LfJ_eval = th.squeeze(batch_dJdx_eval @ fx_eval) # [batch_size, ]
                LgJ_eval = th.squeeze(batch_dJdx_eval @ gx_eval) # [batch_size, ]
                batch_x_square_eval = batch_th_dot(batch_x_eval, batch_x_eval) # [batch_size, ]
                
                batch_h1_eval = (batch_x_square_eval + LfJ_eval - 0.25 * (LgJ_eval**2))    # [batch_size, ]
                
                wrong_number_eval = 0.5*th.sum(th.abs(th.sign(LgJ_eval) + th.sign(batch_lqr_u_eval)))
                sign_loss_eval = 0.01*th.sum(th.relu(LgJ_eval * batch_lqr_u_eval))
                residual_loss_eval = (th.sum(th.abs(batch_h1_eval)))/batch_size_train + th.dot(model(x0), model(x0))
                loss_eval += sign_loss_eval + residual_loss_eval
                
            loss_eval /= 16
            if loss_eval < best_val_loss:
                best_val_loss = loss_eval
                # th.save(model.state_dict(), './CTOptimalControlExperiment/Linear/model/best_model.pth')
                
            print(f'epoch: {i}, loss: {loss.item()}, loss_eval: {loss_eval}')
            # th.save(model.state_dict(), './CTOptimalControlExperiment/Linear/model/last_model.pth')
                
