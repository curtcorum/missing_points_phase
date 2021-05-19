%clean up fid phase difference

close all;

N=32; % nimber of fid pts to use
N_idxs = 1:N;

fid_B = testobj.fid_nav(1:N, 9, 9);
fid_A = testobj.fid_nav(1:N, 9, 1);

%raw difference
fid_phase_A_sub_B = angle( fid_B) - angle( fid_A);
%figure; plot( fid_phase_A_sub_B);

% derivative of difference
%fid_phase_A_sub_B_diff = diff( fid_phase_A_sub_B, 1)
fid_phase_A_sub_B_diff = diff( fid_phase_A_sub_B, 1) + diff( circshift(fid_phase_A_sub_B, -1), 1); % centered derivative
fid_phase_A_sub_B_diff = fid_phase_A_sub_B_diff/2; % centered derivative
%figure; plot( fid_phase_A_sub_B_diff, 'k. ');

% remainder of derivative mod pi/2
fid_phase_A_sub_B_diff = rem( fid_phase_A_sub_B_diff, pi/2);
%figure; plot( fid_phase_A_sub_B_diff_rem);

% remove ouliers
[fid_phase_A_sub_B_diff, removed_idxs] = rmoutliers( fid_phase_A_sub_B_diff);
%removed_idxs = removed_idxs
N_removed = sum( removed_idxs);

% mean is slope of phase, radians/sample, looks good!
b_phase = mean( fid_phase_A_sub_B_diff);
b_phase_std = std( fid_phase_A_sub_B_diff);

N_removed = N_removed   % debug display
b_phase_mean = b_phase  % debug display
b_phase_std = b_phase_std   % debug display

% comparison plots
figure; hold on;
plot( N_idxs, fid_phase_A_sub_B, 'g. ');
plot( N_idxs(removed_idxs), fid_phase_A_sub_B(removed_idxs), 'r. ');
plot( N_idxs, b_phase*(1:N), 'b -');
legend( 'boxon');
if N_removed > 0
legend( 'raw', 'removed', 'fit', 'Location', 'northeastoutside');
else
    legend( 'raw', 'fit', 'Location', 'northeastoutside');
end
%legend( 'raw', 'fit');
hold off;