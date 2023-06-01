
"""
Lowpass filtering the signals : doing 7 for 7/a for now

savedirname = "/net/scratch/lschulz/fluxdata_midwithnee/"*"fluxdata_raw.jld2"
wholedata = SharedArray{Float32}(load(savedirname)["data"])
wholedata_filtered = zeros(Float32,5114,18,16,3)
for i in 1:18
    for j in 1:16
        wholedata_filtered[:,i,j,1] = lowpass_filter(wholedata[:,i,j], 3)
        wholedata_filtered[:,i,j,2] = lowpass_filter(wholedata[:,i,j], 4)
        wholedata_filtered[:,i,j,3] = lowpass_filter(wholedata[:,i,j], 6)

    end
end

jldsave("/net/scratch/lschulz/fluxdata_midwithnee/fluxdata_filtered.jld2",
       data_3a=wholedata_filtered[:,:,:,1],
       data_4a=wholedata_filtered[:,:,:,2],
       data_6a=wholedata_filtered[:,:,:,3])
"""

function lowpass_filter(signal::Vector{T}, cutoff_frequency) where T<:Real
    # Remove the mean of the signal
    signal_mean = mean(signal)
    signal .-= signal_mean

    # Compute the Fourier transform of the signal
    fourier = fft(signal) |> fftshift

    # Compute the frequencies corresponding to the Fourier coefficients
    Ts = 1 / 365.25
    t = 0:Ts:(N-1)*Ts
    freqs = fftfreq(length(t), 1.0/Ts) |> fftshift

    # Find the indices of the Fourier coefficients corresponding to frequencies above the cutoff
    inds = findall(x -> abs(x) > cutoff_frequency, freqs)

    # Zero out the Fourier coefficients corresponding to frequencies above the cutoff
    fourier[inds] .= 0

    # Compute the inverse Fourier transform of the modified Fourier transform
    backtransform = real.(ifft(fourier |> fftshift)) .+ signal_mean

    return real(backtransform)
end
