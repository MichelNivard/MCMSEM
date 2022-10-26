.jit_funcs <- list(
  t4crossprod="
def fn(data, idx2, idx3, idx4):
    data = data - torch.mean(data, dim=0, keepdim=True)
    ncol_out = len(idx2)
    ncol_out += len(idx3) if len(idx3) > 1 else 0
    ncol_out += len(idx4) if len(idx4) > 1 else 0
    out = torch.zeros((data.shape[0], ncol_out), dtype=data.dtype, device=data.device)
    for i, x in enumerate(data):
        x2 = torch.outer(x, x)
        o = [torch.flatten(x2)[idx2]]
        if (len(idx3) > 1) or (len(idx4) > 1):
            x3 = torch.kron(x2, x)
            if len(idx3) > 1:
                o += [torch.flatten(x3)[idx3]]
        else:
            x3 = torch.tensor(0)  # This will never impact anything but is required for compilation
        if len(idx4) > 1:
            x4 = torch.kron(x3, x)
            o += [torch.flatten(x4)[idx4]]
        out[i, :] = torch.cat(o)
    return out
",
  slownecker="
def fn(x, y, kronrow):
    out = torch.zeros((x.shape[0], int(torch.sum(kronrow))), device=x.device)
    firstprod = torch.kron(y, y)  # In case of memory issues: drop this line
    ncol = 0
    for j in range(int(y.shape[1]**3)):
        if kronrow[j]:
            idx0 = j % y.shape[1]
            idx1 = int(j / y.shape[1]) % y.shape[1]
            idx2 = int(j / y.shape[1]**2) % y.shape[1]
            kroncol = torch.kron(firstprod[:, idx2+idx1*y.shape[1]], y[:, idx0]) # In case of memory issues replace with: kroncol = torch.kron(torch.kron(y[:, idx0], y[:, idx1]), y[:, idx2])
            for k in range(out.shape[0]):
                out[k, ncol] = torch.matmul(x[k, :], kroncol)
            ncol += 1
    return out
",
  slowernecker="
def fn(x, y, kronrow):
    out = torch.zeros((x.shape[0], int(torch.sum(kronrow))), device=x.device)
    ncol = 0
    for j in range(int(y.shape[1]**3)):
        if kronrow[j]:
            idx0 = j % y.shape[1]
            idx1 = int(j / y.shape[1]) % y.shape[1]
            idx2 = int(j / y.shape[1]**2) % y.shape[1]
            kroncol = torch.kron(torch.kron(y[:, idx0], y[:, idx1]), y[:, idx2])
            for k in range(out.shape[0]):
                out[k, ncol] = torch.matmul(x[k, :], kroncol)
            ncol += 1
    return out
",
  covlowmem="
# Currently not used, this is a significantly slower but more memory-efficient covariance function
def fn(data):
    data = data - torch.mean(data, dim=1, keepdim=True)
    out = torch.zeros((data.shape[0], data.shape[0]), dtype=data.dtype, device=data.device)
    for i, x in enumerate(data):
        cov = torch.matmul(data[0:i, :], x)
        out[i, 0:i] = cov
        out[0:i, i] = cov
    return out / (data.shape[1] - 1)
",
  getpredmatrices="
# Currently not used, as from my testing it is actually slower compared to the current implementation
# Maybe useful later?
def slownecker(x, y, kronrow):
    out = torch.zeros((x.shape[0], int(torch.sum(kronrow))), device=x.device)
    ncol = 0
    for j in range(int(y.shape[1] ** 3)):
        if kronrow[j] == 1.0:
            idx0 = j % y.shape[1]
            idx1 = int(j / y.shape[1]) % y.shape[1]
            idx2 = int(j / y.shape[1] ** 2) % y.shape[1]
            kroncol = torch.kron(torch.kron(y[:, idx0], y[:, idx1]), y[:, idx2])
            for k in range(out.shape[0]):
                out[k, ncol] = torch.matmul(x[k, :], kroncol)
            ncol += 1
    return out


def fn(par_list: Dict[str, Tensor], torch_masks: Dict[str, Tensor], torch_maps: Dict[str, Tensor], base_matrices: Dict[str, Tensor], use_skewness: List[bool], use_kurtosis: List[bool], diag_s: List[bool],
                           low_memory: List[int]):
    A = torch.add(base_matrices['A'], torch.sum(torch.mul(torch_maps['A'], par_list['A']), dim=2))
    Fm = torch.add(base_matrices['Fm'], torch.sum(torch.mul(torch_maps['Fm'], par_list['Fm']), dim=2))
    S = torch.add(base_matrices['S'], torch.sum(torch.mul(torch_maps['S'], par_list['S']), dim=2))
    if use_skewness:
        Sk = torch.add(base_matrices['Sk'], torch.sum(torch.mul(torch_maps['Sk'], par_list['Sk']), dim=2))
    else:
        Sk = torch.tensor(0)
    if use_kurtosis:
        if diag_s:
            sqrts = torch.sign(S) * torch.sqrt(torch.abs(S))
            skron = torch.kron(torch.diag(sqrts), torch.kron(torch.diag(sqrts), torch.diag(sqrts)))
            K = torch.add(torch.mul(torch.mul(torch.mul(torch.matmul(sqrts, base_matrices['K']), skron), base_matrices['K2']), torch_masks['K']), torch.sum(torch.mul(torch_maps['K'], par_list['K']), dim=2))
        else:
            sqrts = torch.sign(S) * torch.sqrt(torch.abs(S))
            K = torch.add(torch.mul(torch.mul(torch.matmul(torch.matmul(sqrts, base_matrices['K']), torch.kron(sqrts, torch.kron(sqrts, sqrts))), base_matrices['K2']), torch_masks['K']), torch.sum(torch.mul(torch_maps['K'], par_list['K']), dim=2))
    else:
        K = torch.tensor(0)
    M2 = torch.matmul(torch.matmul(torch.matmul(torch.matmul(Fm, torch.inverse(base_matrices['diag_n_p'] - A)), S), torch.transpose(torch.inverse(base_matrices['diag_n_p'] - A), 0, 1)), torch.transpose(Fm, 0, 1))
    if use_skewness:
        M3 = torch.matmul(torch.matmul(torch.matmul(torch.matmul(Fm, torch.inverse(base_matrices['diag_n_p'] - A)), Sk), torch.kron(torch.transpose(torch.inverse(base_matrices['diag_n_p'] - A), 0, 1), torch.transpose(torch.inverse(base_matrices['diag_n_p'] - A), 0, 1))), torch.kron(torch.transpose(Fm, 0, 1), torch.transpose(Fm, 0, 1)))
    else:
        M3 = torch.tensor(0)
    if use_kurtosis:
        if low_memory[0] > 1:
            fmkronrow =  torch.kron(torch.kron(torch.sum(Fm, dim=0), torch.sum(Fm, dim=0)), torch.sum(Fm, dim=0))
            M4 = slownecker(torch.matmul(torch.matmul(Fm, torch.inverse(base_matrices['diag_n_p'] - A)), K), torch.transpose(torch.inverse(base_matrices['diag_n_p'] - A), 0, 1), fmkronrow)
        else:
            M4 = torch.matmul(torch.matmul(torch.matmul(Fm, torch.inverse(base_matrices['diag_n_p'] - A)), K), torch.matmul(torch.kron(torch.kron(torch.transpose(torch.inverse(base_matrices['diag_n_p'] - A), 0, 1), torch.transpose(torch.inverse(base_matrices['diag_n_p'] - A), 0, 1)), torch.transpose(torch.inverse(base_matrices['diag_n_p'] - A), 0, 1)), torch.kron(torch.kron(torch.transpose(Fm, 0, 1), torch.transpose(Fm, 0, 1)), torch.transpose(Fm, 0, 1))))
    else:
        M4 = torch.tensor(0)

    if use_kurtosis and use_skewness:
        return dict(M2=M2, M3=M3, M4=M4)
    elif use_skewness:
        return dict(M2=M2, M3=M3)
    elif use_kurtosis:
        return dict(M2=M2, M4=M4)
    else:
        return dict(M2=M2)
"
)