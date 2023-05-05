function [gradient] = LargerGradient(im_s, im_t, idx_i, idx_j)

n_i = idx_i(1);
m_i = idx_i(2);
c_i = idx_i(3);

n_j = idx_j(1);
m_j = idx_j(2);
c_j = idx_j(3);

g_s = im_s(n_i,m_i,c_i) - im_s(n_j,m_j,c_j);
g_t = im_t(n_i,m_i,c_i) - im_t(n_j,m_j,c_j);

if abs(g_s) < abs(g_t)
    gradient = g_t;
else
    gradient = g_s;
end

