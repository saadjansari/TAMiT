function w = ProxHard(T, w)

w = w .* (abs(w) < T);