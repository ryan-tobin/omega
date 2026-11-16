//! Work on the atoms split across threads.
//!
//! The items, one per atom, are cut into runs of neighboring items, one
//! run per thread, and each thread changes only the items of its own run.
//! No two threads touch the same item, so nothing needs a lock, and each
//! item gets the same result whatever the number of threads. The threads
//! come from the standard library's std::thread::scope, which waits for
//! every thread it started before it returns.

use std::thread;

/// Calls `work(i, &mut items[i])` once for every index i. The items are
/// split into at most `threads` runs of neighboring items, each handled by
/// its own thread. The calling thread takes the first run itself, so one
/// thread starts no other. `threads` must be at least 1.
pub fn for_each<T, F>(items: &mut [T], threads: usize, work: F)
where
    T: Send,
    F: Fn(usize, &mut T) + Sync,
{
    assert!(threads >= 1, "need at least one thread");
    let size = items.len().div_ceil(threads).max(1);
    let work = &work;
    thread::scope(|scope| {
        let mut runs = items.chunks_mut(size).enumerate();
        let first = runs.next();
        for (k, run) in runs {
            scope.spawn(move || each(k * size, run, work));
        }
        if let Some((_, run)) = first {
            each(0, run, work);
        }
    });
}

/// Calls `work` on each item of `run`, whose first item has index `start`.
fn each<T, F: Fn(usize, &mut T)>(start: usize, run: &mut [T], work: &F) {
    for (k, item) in run.iter_mut().enumerate() {
        work(start + k, item);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    #[test]
    fn every_item_is_visited_once_with_its_own_index() {
        // Each item records the index it was given and how many times it
        // was visited. Exact: (i, 1) for every i, for any number of items
        // and threads, including no items and more threads than items.
        for n in [0, 1, 5, 64, 1000] {
            for threads in [1, 2, 3, 4, 7, 16] {
                let mut items = vec![(usize::MAX, 0); n];
                for_each(&mut items, threads, |i, item| *item = (i, item.1 + 1));
                for (i, &item) in items.iter().enumerate() {
                    assert_eq!(item, (i, 1), "{n} items, {threads} threads");
                }
            }
        }
    }

    #[test]
    fn work_runs_on_as_many_threads_as_asked() {
        // 100 items on 4 threads: runs of 25, each on its own thread, the
        // first on the calling one.
        let mut items = vec![None; 100];
        for_each(&mut items, 4, |_, item| {
            *item = Some(thread::current().id())
        });
        let ids: HashSet<_> = items.iter().map(|id| id.unwrap()).collect();
        assert_eq!(ids.len(), 4);
        assert!(ids.contains(&thread::current().id()));
    }

    #[test]
    #[should_panic(expected = "need at least one thread")]
    fn zero_threads_are_refused() {
        for_each(&mut [0.0; 3], 0, |_, x| *x += 1.0);
    }
}