class Video:
    def __init__(self, url: str, caption: str = "") -> None:
        self.url = url
        self.caption = caption

    def as_markdown(self) -> str:
        if not self.url or self.url == "":
            return "![]()"
        return f"![{self.caption}]({self.format_url()})"

    def format_url(self) -> str:
        if self.url.endswith(".mp4"):
            return self.url

        return self.url + ".mp4"
