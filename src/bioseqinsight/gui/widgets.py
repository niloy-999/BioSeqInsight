"""Small reusable Tkinter widgets.

Keeping these separate means the view modules contain layout and delegation
only. No widget in this package performs a biological calculation or an HTTP
request; they render what the core API returns.
"""

from __future__ import annotations

import queue
import threading
import tkinter as tk
from collections.abc import Callable
from tkinter import ttk
from typing import Any

MONO = ("Menlo", 11) if tk.TkVersion >= 8.6 else ("Courier", 11)

MAPPING_COLOURS = {
    "M4": "#1a7f37",
    "M3": "#0969da",
    "M2": "#9a6700",
    "M1": "#cf222e",
    "M0": "#6e7781",
}


class LabelledText(ttk.Frame):
    """A labelled, scrollable text area used for input and output panes."""

    def __init__(self, master, label: str, height: int = 8, readonly: bool = False, **kwargs):
        super().__init__(master, **kwargs)
        self.columnconfigure(0, weight=1)
        self.rowconfigure(1, weight=1)
        ttk.Label(self, text=label).grid(row=0, column=0, sticky="w", pady=(0, 4))

        container = ttk.Frame(self)
        container.grid(row=1, column=0, sticky="nsew")
        container.columnconfigure(0, weight=1)
        container.rowconfigure(0, weight=1)

        self.text = tk.Text(container, height=height, wrap="none", font=MONO, undo=True)
        self.text.grid(row=0, column=0, sticky="nsew")
        y_scroll = ttk.Scrollbar(container, orient="vertical", command=self.text.yview)
        y_scroll.grid(row=0, column=1, sticky="ns")
        x_scroll = ttk.Scrollbar(container, orient="horizontal", command=self.text.xview)
        x_scroll.grid(row=1, column=0, sticky="ew")
        self.text.configure(yscrollcommand=y_scroll.set, xscrollcommand=x_scroll.set)
        self.readonly = readonly
        if readonly:
            self.text.configure(state="disabled")

    def get(self) -> str:
        return self.text.get("1.0", "end-1c")

    def set(self, value: str) -> None:
        if self.readonly:
            self.text.configure(state="normal")
        self.text.delete("1.0", "end")
        self.text.insert("1.0", value)
        if self.readonly:
            self.text.configure(state="disabled")

    def append(self, value: str) -> None:
        if self.readonly:
            self.text.configure(state="normal")
        self.text.insert("end", value)
        self.text.see("end")
        if self.readonly:
            self.text.configure(state="disabled")

    def clear(self) -> None:
        self.set("")


class StatusBar(ttk.Frame):
    """Status line with an indeterminate progress bar for background work."""

    def __init__(self, master, **kwargs):
        super().__init__(master, **kwargs)
        self.columnconfigure(0, weight=1)
        self._message = tk.StringVar(value="Ready.")
        ttk.Label(self, textvariable=self._message, anchor="w").grid(row=0, column=0, sticky="ew")
        self.progress = ttk.Progressbar(self, mode="indeterminate", length=140)
        self.progress.grid(row=0, column=1, sticky="e", padx=(8, 0))
        self.progress.grid_remove()

    def set(self, message: str) -> None:
        self._message.set(message)

    def busy(self, message: str = "Working...") -> None:
        self.set(message)
        self.progress.grid()
        self.progress.start(12)

    def idle(self, message: str = "Ready.") -> None:
        self.progress.stop()
        self.progress.grid_remove()
        self.set(message)


class BackgroundRunner:
    """Run a callable off the Tk event loop and deliver the result safely.

    Network calls must never run on the main thread or the whole window
    freezes while a structural service times out. Tk widgets must never be
    touched from a worker thread either, so results come back through a queue
    that the main loop polls.
    """

    def __init__(self, widget: tk.Misc, poll_ms: int = 80):
        self.widget = widget
        self.poll_ms = poll_ms
        self._queue: queue.Queue[tuple[Callable, tuple]] = queue.Queue()
        self._polling = False

    def submit(
        self,
        work: Callable[[], Any],
        on_success: Callable[[Any], None],
        on_error: Callable[[Exception], None] | None = None,
    ) -> threading.Thread:
        def target() -> None:
            try:
                result = work()
            except Exception as exc:
                self._queue.put((on_error or _reraise, (exc,)))
            else:
                self._queue.put((on_success, (result,)))

        thread = threading.Thread(target=target, daemon=True)
        thread.start()
        self._ensure_polling()
        return thread

    def _ensure_polling(self) -> None:
        if self._polling:
            return
        self._polling = True
        self._poll()

    def _poll(self) -> None:
        drained = False
        while True:
            try:
                callback, args = self._queue.get_nowait()
            except queue.Empty:
                break
            drained = True
            try:
                callback(*args)
            except Exception:  # pragma: no cover - callback bugs must not kill the loop
                import traceback

                traceback.print_exc()
        self._polling = bool(not drained or True)
        self.widget.after(self.poll_ms, self._poll)


def _reraise(exc: Exception) -> None:  # pragma: no cover - default error path
    raise exc


class MappingBadge(ttk.Label):
    """Coloured badge showing an M0-M4 mapping level."""

    def __init__(self, master, **kwargs):
        super().__init__(master, text="", **kwargs)

    def show(self, level_value: str, label: str) -> None:
        self.configure(
            text=f"{level_value} - {label}",
            foreground=MAPPING_COLOURS.get(level_value, "#6e7781"),
        )


def bind_copy(text_widget: tk.Text) -> None:
    """Allow select-all and copy in read-only output panes."""

    def select_all(event):
        event.widget.tag_add("sel", "1.0", "end")
        return "break"

    text_widget.bind("<Control-a>", select_all)
    text_widget.bind("<Command-a>", select_all)


__all__ = [
    "MAPPING_COLOURS",
    "BackgroundRunner",
    "LabelledText",
    "MappingBadge",
    "StatusBar",
    "bind_copy",
]
