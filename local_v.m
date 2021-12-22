function v_local = local_v(vp,c_eff,c_plat,v0)

% this function calculates the local velocity vp (in the main program: 
% "v(pp)") of particle pp based on the effective concentration v_eff "how
% fast is the particle actually?", the plateau concentration (if we want a
% plateau in the profile) and a certain scaling concentration v0.
% only used for Chemokinesis.

%%note that the plateau was rejected during the process.

if c_eff > c_plat
    v_local = v0 + vp;
else
    v_local = v0 + (c_eff-c_plat)/c_plat.*v0 + vp;
end

end
