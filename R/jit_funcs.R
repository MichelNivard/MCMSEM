.jit_funcs <- list(
  t4crossprod="
def fn(data, idx2, idx3, idx4):
    data = data - torch.mean(data, dim=0, keepdim=True)
    out = torch.zeros((len(idx2)+len(idx3)+len(idx4), data.shape[0]), dtype=data.dtype, device=data.device)
    for i, x in enumerate(data):
        x2 = torch.outer(x, x)
        x3 = torch.kron(x2, x)
        x4 = torch.kron(x3, x)
        out[:, i] = torch.cat([torch.flatten(x2)[idx2], torch.flatten(x3)[idx3], torch.flatten(x4)[idx4]])
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
"
)