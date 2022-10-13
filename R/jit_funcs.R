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
"
)